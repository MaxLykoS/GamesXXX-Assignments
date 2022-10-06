#ifdef GL_ES
precision highp float;
#endif

uniform vec3 uLightDir;
uniform vec3 uCameraPos;
uniform vec3 uLightRadiance;
uniform sampler2D uGDiffuse;
uniform sampler2D uGDepth;
uniform sampler2D uGNormalWorld;
uniform sampler2D uGShadow;
uniform sampler2D uGPosWorld;

varying mat4 vWorldToScreen;
varying highp vec4 vPosWorld;

#define M_PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307
#define INV_PI 0.31830988618
#define INV_TWO_PI 0.15915494309

float Rand1(inout float p) {
  p = fract(p * .1031);
  p *= p + 33.33;
  p *= p + p;
  return fract(p);
}

vec2 Rand2(inout float p) {
  return vec2(Rand1(p), Rand1(p));
}

float InitRand(vec2 uv) {
	vec3 p3  = fract(vec3(uv.xyx) * .1031);
  p3 += dot(p3, p3.yzx + 33.33);
  return fract((p3.x + p3.y) * p3.z);
}

vec3 SampleHemisphereUniform(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = uv.x;
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(1.0 - z*z);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = INV_TWO_PI;
  return dir;
}

vec3 SampleHemisphereCos(inout float s, out float pdf) {
  vec2 uv = Rand2(s);
  float z = sqrt(1.0 - uv.x);
  float phi = uv.y * TWO_PI;
  float sinTheta = sqrt(uv.x);
  vec3 dir = vec3(sinTheta * cos(phi), sinTheta * sin(phi), z);
  pdf = z * INV_PI;
  return dir;
}

void LocalBasis(vec3 n, out vec3 b1, out vec3 b2) {
  float sign_ = sign(n.z);
  if (n.z == 0.0) {
    sign_ = 1.0;
  }
  float a = -1.0 / (sign_ + n.z);
  float b = n.x * n.y * a;
  b1 = vec3(1.0 + sign_ * n.x * n.x * a, sign_ * b, -sign_ * n.x);
  b2 = vec3(b, sign_ + n.y * n.y * a, -n.y);
}

vec4 Project(vec4 a) {
  return a / a.w;
}

float GetDepth(vec3 posWorld) {
  float depth = (vWorldToScreen * vec4(posWorld, 1.0)).w;
  return depth;
}

/*
 * Transform point from world space to screen space([0, 1] x [0, 1])
 *
 */
vec2 GetScreenCoordinate(vec3 posWorld) {
  vec2 uv = Project(vWorldToScreen * vec4(posWorld, 1.0)).xy * 0.5 + 0.5;
  return uv;
}

float GetGBufferDepth(vec2 uv) {
  float depth = texture2D(uGDepth, uv).x;
  if (depth < 1e-2) {
    depth = 1000.0;
  }
  return depth;
}

vec3 GetGBufferNormalWorld(vec2 uv) {
  vec3 normal = texture2D(uGNormalWorld, uv).xyz;
  return normal;
}

vec3 GetGBufferPosWorld(vec2 uv) {
  vec3 posWorld = texture2D(uGPosWorld, uv).xyz;
  return posWorld;
}

float GetGBufferuShadow(vec2 uv) {
  float visibility = texture2D(uGShadow, uv).x;
  return visibility;
}

vec3 GetGBufferDiffuse(vec2 uv) {
  vec3 diffuse = texture2D(uGDiffuse, uv).xyz;
  diffuse = pow(diffuse, vec3(2.2));
  return diffuse;
}

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;

    float nom   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = M_PI * denom * denom;

    return nom / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 EvalDiffuse(vec3 V, vec3 L, vec2 screenUV)
{
  vec3 Lo = vec3(0.0);

  vec3 albedo = GetGBufferDiffuse(screenUV);
  vec3 worldPos = GetGBufferPosWorld(screenUV);
  vec3 N = GetGBufferNormalWorld(screenUV);
  float depth = GetGBufferDepth(screenUV);

  // cook-torrance BRDF
  float roughness = 0.3;
  float metallic = 0.0;
  vec3 H = normalize(V + L);

  vec3 F0 = vec3(0.04); 
  F0 = mix(F0, albedo, metallic);

  float NDF = DistributionGGX(N, H, roughness);        
  float G = GeometrySmith(N, V, L, roughness);      
  vec3 F = fresnelSchlick(max(dot(H, V), 0.0), F0);       

  vec3 kS = F;
  vec3 kD = vec3(1.0) - kS;
  kD *= 1.0 - metallic;     

  vec3 nominator = NDF * G * F;
  float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.001; 
  vec3 specular = nominator / denominator;

  // add to outgoing radiance Lo
  float NdotL = max(dot(N, L), 0.0);                
  Lo += (kD * albedo / M_PI + specular) * uLightRadiance * NdotL;

  return Lo;
}
vec3 EvalDirectionalLight(vec2 uv) 
{
  vec3 posW = GetGBufferPosWorld(uv);
  vec3 wi = normalize(uLightDir);
  vec3 wo = normalize(uCameraPos - posW);
  vec3 bsdf = EvalDiffuse(wi, wo, uv);
  return uLightRadiance * bsdf * GetGBufferuShadow(uv);
}

bool outScreen(vec3 pos)
{
  vec2 uv = GetScreenCoordinate(pos);
  return any(bvec4(lessThan(uv, vec2(0.0)), greaterThan(uv, vec2(1.0))));
}
bool RayMarch(vec3 ori, vec3 dir, out vec3 hitPos) 
{
  vec3 cur = ori;
  dir = normalize(dir);
  for(int i = 0; i < 100; ++i)
  {
    if(outScreen(cur))
      break;
    vec2 uv = GetScreenCoordinate(cur);
    float hitDepth = GetGBufferDepth(uv);
    float rayDepth = GetDepth(cur);
    if(rayDepth - 0.01 > hitDepth)
    {
      hitPos = cur;
      return true;
    }
    else
      cur = cur + dir * 0.5;
  }
  return false;
}
#define INIT_STEP 0.8
#define MAX_STEPS 20
#define EPS 1e-2
#define THRES 0.1
bool atFront(vec3 pos){
  return GetDepth(pos) < GetGBufferDepth(GetScreenCoordinate(pos));
}
bool hasInter(vec3 pos, vec3 dir, out vec3 hitPos){
  float d1 = GetGBufferDepth(GetScreenCoordinate(pos)) - GetDepth(pos) + EPS;
  float d2 = GetDepth(pos + dir) - GetGBufferDepth(GetScreenCoordinate(pos + dir)) + EPS;
  if(d1 < THRES && d2 < THRES){
    hitPos = pos + dir * d1 / (d1 + d2);
    return true;
  }  
  return false;
}

bool RayMarchh(vec3 ori, vec3 dir, out vec3 hitPos) {
  bool intersect = false, firstinter = false;
  float st = INIT_STEP;
  vec3 current = ori;
  for (int i = 0;i < MAX_STEPS;i++){
    if(outScreen(current)){
      break;
    }
    else if(atFront(current + dir * st)){
      current += dir * st;
    }else{
      firstinter = true;
      if(st < EPS){
        if(hasInter(current, dir * st * 2.0, hitPos)){
          intersect = true;
        }
        break;
      }
    }
    if(firstinter)
      st *= 0.5;
  }
  return intersect;
}

#define SAMPLE_NUM 4
vec3 EvalIndirectLight(vec3 pos)
{
  float pdf, seed = dot(pos, vec3(100.0));
  vec3 Li = vec3(0.0), dir, hitPos;
  vec3 normal = GetGBufferNormalWorld(GetScreenCoordinate(pos)), b1, b2;
  LocalBasis(normal, b1, b2);
  mat3 TBN = mat3(b1, b2, normal);
  for(int i = 0; i < SAMPLE_NUM;i++)
  {
    dir = normalize(TBN * SampleHemisphereCos(seed, pdf));
    if(RayMarchh(pos, dir, hitPos))
    {
      vec3 wo = normalize(uCameraPos - pos);
      vec3 L = EvalDiffuse(dir, wo, GetScreenCoordinate(pos)) / pdf;
      wo = normalize(uCameraPos - hitPos);
      vec3 wi = normalize(uLightDir);
      L *= EvalDiffuse(wi, wo, GetScreenCoordinate(hitPos)) * EvalDirectionalLight(GetScreenCoordinate(hitPos));
      Li += L;
    }
  }
  return Li / float(SAMPLE_NUM);
}

void main() 
{
  float s = InitRand(gl_FragCoord.xy);
  vec3 Lo = vec3(0.0,0.0,0.0);

  vec2 screenUV = GetScreenCoordinate(vPosWorld.xyz);
  vec3 V = normalize(uCameraPos - GetGBufferPosWorld(screenUV));
  vec3 L = normalize(uLightDir);

  Lo = EvalDiffuse(V, L, screenUV);
  Lo *= GetGBufferuShadow(screenUV);
  Lo = Lo + EvalIndirectLight(vPosWorld.xyz)*10.0;

  vec3 color = pow(clamp(Lo, vec3(0.0), vec3(1.0)), vec3(1.0 / 2.2));
  gl_FragColor = vec4(vec3(color.rgb), 1.0);
}
