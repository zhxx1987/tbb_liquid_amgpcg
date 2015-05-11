/*
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#define STRINGIFY(A) #A

// vertex shader
const char *vertexShader = STRINGIFY(
                               uniform float pointRadius;  // point size in world space
                               uniform float pointScale;   // scale to calculate size in pixels
                               uniform float densityScale;
                               uniform float densityOffset;
							   varying vec3 v;
							   varying float r;
                               void main()
{
    // calculate window-space point size
    vec3 posEye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
	v = posEye;
	r = pointRadius;
    float dist = length(posEye);
    gl_PointSize = pointRadius * (pointScale / dist);

    gl_TexCoord[0] = gl_MultiTexCoord0;
	vec3 pos = gl_Vertex.xyz;
	if (pos.z>0.25)
	{ 
		pos = vec3(0,0,0);
	}
    gl_Position = gl_ModelViewProjectionMatrix * vec4(pos, 1.0);

    gl_FrontColor = vec4(0.1,0.2,0.8,1.0);
}
                           );

// pixel shader for rendering points as shaded spheres
const char *spherePixelShader = STRINGIFY(
	varying vec3 v;
	varying float r;
                                    void main()
{
	
    const vec3 lightDir = vec3(0.577, 0.577, 0.577);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(N.xy, N.xy);

    if (mag > 1.0) discard;   // kill pixels outside circle

    N.z = sqrt(1.0-mag);

	vec3 pv = vec3(v.x,v.y,v.z-N.z*r);

	vec3 E = normalize(-pv);
	vec3 R = normalize(-reflect(lightDir,N)); 
	vec4 Ispec = 1.0 * pow(max(dot(R,E),0.0),1.3);
	Ispec = clamp(Ispec, 0.0, 1.0); 
    // calculate lighting
    float diffuse = max(0.0, dot(lightDir, N));
	vec4 color = vec4(0.1,0.2,0.8,1.0);
    gl_FragColor = color*diffuse + color*Ispec;
}
                                );
