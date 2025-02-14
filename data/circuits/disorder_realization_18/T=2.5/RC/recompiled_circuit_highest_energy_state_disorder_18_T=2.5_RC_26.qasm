OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9189605) q[0];
sx q[0];
rz(-0.55813342) q[0];
sx q[0];
rz(-2.4827935) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(0.08858362) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9889116) q[0];
sx q[0];
rz(-1.7636239) q[0];
sx q[0];
rz(1.1367984) q[0];
x q[1];
rz(-2.2031783) q[2];
sx q[2];
rz(-0.66667306) q[2];
sx q[2];
rz(-3.0423903) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7329053) q[1];
sx q[1];
rz(-0.63893203) q[1];
sx q[1];
rz(-1.8270135) q[1];
rz(-pi) q[2];
rz(-2.5837059) q[3];
sx q[3];
rz(-2.2672664) q[3];
sx q[3];
rz(-0.59244746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6170071) q[2];
sx q[2];
rz(-1.1716537) q[2];
sx q[2];
rz(-0.48801547) q[2];
rz(0.46767849) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(-2.974143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92983285) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(-1.8875341) q[0];
rz(2.5414741) q[1];
sx q[1];
rz(-1.80872) q[1];
sx q[1];
rz(2.7235203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.583823) q[0];
sx q[0];
rz(-2.1033786) q[0];
sx q[0];
rz(1.363165) q[0];
rz(-pi) q[1];
rz(0.64855021) q[2];
sx q[2];
rz(-1.4709215) q[2];
sx q[2];
rz(1.5431736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4198735) q[1];
sx q[1];
rz(-1.5054107) q[1];
sx q[1];
rz(-1.2544778) q[1];
rz(-pi) q[2];
rz(-2.3206059) q[3];
sx q[3];
rz(-2.0625522) q[3];
sx q[3];
rz(3.0420835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40689251) q[2];
sx q[2];
rz(-1.1310581) q[2];
sx q[2];
rz(-3.0231754) q[2];
rz(-2.7393869) q[3];
sx q[3];
rz(-1.4879358) q[3];
sx q[3];
rz(1.8124628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44620946) q[0];
sx q[0];
rz(-2.8948247) q[0];
sx q[0];
rz(-1.1545908) q[0];
rz(2.0788976) q[1];
sx q[1];
rz(-2.1176391) q[1];
sx q[1];
rz(-1.9564995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20511928) q[0];
sx q[0];
rz(-1.8471878) q[0];
sx q[0];
rz(-2.2689681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22151557) q[2];
sx q[2];
rz(-2.0755526) q[2];
sx q[2];
rz(-1.5919478) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.628988) q[1];
sx q[1];
rz(-1.1956685) q[1];
sx q[1];
rz(0.30156044) q[1];
rz(-pi) q[2];
rz(0.55425476) q[3];
sx q[3];
rz(-0.46559428) q[3];
sx q[3];
rz(-1.3265644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9065173) q[2];
sx q[2];
rz(-0.86248988) q[2];
sx q[2];
rz(-1.6974576) q[2];
rz(-2.2200572) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(-0.055189341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7017355) q[0];
sx q[0];
rz(-0.12049645) q[0];
sx q[0];
rz(-3.0635656) q[0];
rz(1.1174508) q[1];
sx q[1];
rz(-1.2470587) q[1];
sx q[1];
rz(1.6695581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271281) q[0];
sx q[0];
rz(-2.3660866) q[0];
sx q[0];
rz(1.7239611) q[0];
rz(2.1536768) q[2];
sx q[2];
rz(-1.1088409) q[2];
sx q[2];
rz(-0.5505901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8824655) q[1];
sx q[1];
rz(-1.5423911) q[1];
sx q[1];
rz(-1.9309429) q[1];
rz(1.1999287) q[3];
sx q[3];
rz(-2.2590593) q[3];
sx q[3];
rz(2.1394068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3052519) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(0.15023896) q[2];
rz(-2.5495106) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(1.1881812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2109461) q[0];
sx q[0];
rz(-0.35363126) q[0];
sx q[0];
rz(2.3760702) q[0];
rz(-2.3402479) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(1.9498922) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0203637) q[0];
sx q[0];
rz(-2.1018711) q[0];
sx q[0];
rz(-1.2174106) q[0];
x q[1];
rz(1.0982047) q[2];
sx q[2];
rz(-2.049438) q[2];
sx q[2];
rz(0.20028534) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5687823) q[1];
sx q[1];
rz(-2.2523419) q[1];
sx q[1];
rz(3.0589025) q[1];
rz(-pi) q[2];
rz(-1.839961) q[3];
sx q[3];
rz(-1.3084305) q[3];
sx q[3];
rz(-0.056726168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8754862) q[2];
sx q[2];
rz(-0.85662872) q[2];
sx q[2];
rz(2.6066656) q[2];
rz(0.62355012) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(2.258544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542338) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(2.2075388) q[0];
rz(1.2454698) q[1];
sx q[1];
rz(-1.5555236) q[1];
sx q[1];
rz(2.5159871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4564292) q[0];
sx q[0];
rz(-2.3722088) q[0];
sx q[0];
rz(-1.3145694) q[0];
rz(1.6883739) q[2];
sx q[2];
rz(-2.46799) q[2];
sx q[2];
rz(-0.51076159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4491419) q[1];
sx q[1];
rz(-2.3896273) q[1];
sx q[1];
rz(2.7099931) q[1];
rz(-pi) q[2];
rz(-0.91208338) q[3];
sx q[3];
rz(-0.28277031) q[3];
sx q[3];
rz(1.2706437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.049456747) q[2];
sx q[2];
rz(-2.226604) q[2];
sx q[2];
rz(0.12506872) q[2];
rz(-0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(0.47621134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8759988) q[0];
sx q[0];
rz(-0.47530526) q[0];
sx q[0];
rz(-1.2566316) q[0];
rz(1.9427293) q[1];
sx q[1];
rz(-1.7190944) q[1];
sx q[1];
rz(1.2748324) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0755646) q[0];
sx q[0];
rz(-0.69113737) q[0];
sx q[0];
rz(0.71623556) q[0];
x q[1];
rz(-0.9918757) q[2];
sx q[2];
rz(-2.5235368) q[2];
sx q[2];
rz(1.2787873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82060888) q[1];
sx q[1];
rz(-2.5638691) q[1];
sx q[1];
rz(-2.251365) q[1];
x q[2];
rz(-1.8883287) q[3];
sx q[3];
rz(-1.4503363) q[3];
sx q[3];
rz(0.61653709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4510497) q[2];
sx q[2];
rz(-2.1284911) q[2];
sx q[2];
rz(-0.96735442) q[2];
rz(1.5830154) q[3];
sx q[3];
rz(-1.5019006) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3485182) q[0];
sx q[0];
rz(-0.59978849) q[0];
sx q[0];
rz(-3.0317958) q[0];
rz(1.173136) q[1];
sx q[1];
rz(-2.1497612) q[1];
sx q[1];
rz(1.0822302) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8380676) q[0];
sx q[0];
rz(-0.81373022) q[0];
sx q[0];
rz(-2.8007048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3668109) q[2];
sx q[2];
rz(-0.49353853) q[2];
sx q[2];
rz(-0.94369027) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6746862) q[1];
sx q[1];
rz(-1.9787346) q[1];
sx q[1];
rz(-0.80927421) q[1];
rz(-0.63858648) q[3];
sx q[3];
rz(-1.2632086) q[3];
sx q[3];
rz(1.5945376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75862306) q[2];
sx q[2];
rz(-2.0489645) q[2];
sx q[2];
rz(0.36337241) q[2];
rz(-2.162497) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30524224) q[0];
sx q[0];
rz(-0.79638052) q[0];
sx q[0];
rz(2.7889732) q[0];
rz(-1.7800219) q[1];
sx q[1];
rz(-1.9510599) q[1];
sx q[1];
rz(-0.1415267) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3985719) q[0];
sx q[0];
rz(-1.6048661) q[0];
sx q[0];
rz(2.9719246) q[0];
rz(-pi) q[1];
rz(-2.5078819) q[2];
sx q[2];
rz(-0.6978718) q[2];
sx q[2];
rz(-0.77264589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80672164) q[1];
sx q[1];
rz(-1.1523243) q[1];
sx q[1];
rz(-1.1711907) q[1];
rz(-pi) q[2];
rz(2.4961595) q[3];
sx q[3];
rz(-1.107599) q[3];
sx q[3];
rz(0.86477676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8286459) q[2];
sx q[2];
rz(-1.4742278) q[2];
sx q[2];
rz(2.055577) q[2];
rz(-0.18276754) q[3];
sx q[3];
rz(-2.1153617) q[3];
sx q[3];
rz(-2.1695547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(-2.4973448) q[0];
rz(3.0746025) q[1];
sx q[1];
rz(-1.7552152) q[1];
sx q[1];
rz(3.0279874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035126) q[0];
sx q[0];
rz(-1.3467731) q[0];
sx q[0];
rz(-3.1372848) q[0];
rz(-pi) q[1];
x q[1];
rz(3.058931) q[2];
sx q[2];
rz(-1.4689406) q[2];
sx q[2];
rz(1.5993303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1190824) q[1];
sx q[1];
rz(-1.4003229) q[1];
sx q[1];
rz(-1.4524231) q[1];
rz(-pi) q[2];
rz(1.6842277) q[3];
sx q[3];
rz(-0.3598752) q[3];
sx q[3];
rz(0.42478334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1429448) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(-2.3764853) q[2];
rz(0.1782002) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(-2.3044738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2160303) q[0];
sx q[0];
rz(-2.3860274) q[0];
sx q[0];
rz(-0.26269333) q[0];
rz(-2.8021011) q[1];
sx q[1];
rz(-1.9120293) q[1];
sx q[1];
rz(1.0614352) q[1];
rz(-0.65503623) q[2];
sx q[2];
rz(-1.8476386) q[2];
sx q[2];
rz(2.4551433) q[2];
rz(2.8398561) q[3];
sx q[3];
rz(-1.269125) q[3];
sx q[3];
rz(2.453809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
