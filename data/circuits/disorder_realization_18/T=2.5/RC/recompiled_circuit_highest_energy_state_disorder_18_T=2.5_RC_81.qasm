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
rz(0.65879917) q[0];
rz(0.88762033) q[1];
sx q[1];
rz(-0.69488156) q[1];
sx q[1];
rz(0.08858362) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6348931) q[0];
sx q[0];
rz(-1.1453712) q[0];
sx q[0];
rz(-0.21197196) q[0];
rz(-1.0052571) q[2];
sx q[2];
rz(-1.1966248) q[2];
sx q[2];
rz(-2.1924874) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.772208) q[1];
sx q[1];
rz(-1.7225035) q[1];
sx q[1];
rz(-0.94767344) q[1];
x q[2];
rz(-1.0064355) q[3];
sx q[3];
rz(-2.2794323) q[3];
sx q[3];
rz(-0.17911994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6170071) q[2];
sx q[2];
rz(-1.969939) q[2];
sx q[2];
rz(-0.48801547) q[2];
rz(2.6739142) q[3];
sx q[3];
rz(-0.52507639) q[3];
sx q[3];
rz(-0.16744965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92983285) q[0];
sx q[0];
rz(-0.81120315) q[0];
sx q[0];
rz(-1.8875341) q[0];
rz(-2.5414741) q[1];
sx q[1];
rz(-1.80872) q[1];
sx q[1];
rz(-2.7235203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.583823) q[0];
sx q[0];
rz(-1.038214) q[0];
sx q[0];
rz(-1.7784276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64855021) q[2];
sx q[2];
rz(-1.4709215) q[2];
sx q[2];
rz(-1.5431736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1723056) q[1];
sx q[1];
rz(-1.2551771) q[1];
sx q[1];
rz(0.068788485) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3206059) q[3];
sx q[3];
rz(-2.0625522) q[3];
sx q[3];
rz(3.0420835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40689251) q[2];
sx q[2];
rz(-1.1310581) q[2];
sx q[2];
rz(0.1184173) q[2];
rz(0.40220574) q[3];
sx q[3];
rz(-1.6536568) q[3];
sx q[3];
rz(-1.8124628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6953832) q[0];
sx q[0];
rz(-0.24676794) q[0];
sx q[0];
rz(1.9870019) q[0];
rz(-1.062695) q[1];
sx q[1];
rz(-1.0239536) q[1];
sx q[1];
rz(-1.1850932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9364734) q[0];
sx q[0];
rz(-1.8471878) q[0];
sx q[0];
rz(2.2689681) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22151557) q[2];
sx q[2];
rz(-2.0755526) q[2];
sx q[2];
rz(-1.5496448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.628988) q[1];
sx q[1];
rz(-1.9459241) q[1];
sx q[1];
rz(-0.30156044) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55425476) q[3];
sx q[3];
rz(-0.46559428) q[3];
sx q[3];
rz(-1.3265644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2350754) q[2];
sx q[2];
rz(-0.86248988) q[2];
sx q[2];
rz(1.4441351) q[2];
rz(0.92153543) q[3];
sx q[3];
rz(-2.5369365) q[3];
sx q[3];
rz(-0.055189341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4398572) q[0];
sx q[0];
rz(-0.12049645) q[0];
sx q[0];
rz(-3.0635656) q[0];
rz(-1.1174508) q[1];
sx q[1];
rz(-1.8945339) q[1];
sx q[1];
rz(-1.4720346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8142151) q[0];
sx q[0];
rz(-0.80670415) q[0];
sx q[0];
rz(-2.9931158) q[0];
rz(-pi) q[1];
rz(2.3061137) q[2];
sx q[2];
rz(-2.4148921) q[2];
sx q[2];
rz(1.5269321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2591272) q[1];
sx q[1];
rz(-1.5992016) q[1];
sx q[1];
rz(-1.9309429) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7265161) q[3];
sx q[3];
rz(-2.3743694) q[3];
sx q[3];
rz(-1.5900714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3052519) q[2];
sx q[2];
rz(-0.64771104) q[2];
sx q[2];
rz(2.9913537) q[2];
rz(2.5495106) q[3];
sx q[3];
rz(-1.3034857) q[3];
sx q[3];
rz(-1.1881812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2109461) q[0];
sx q[0];
rz(-2.7879614) q[0];
sx q[0];
rz(2.3760702) q[0];
rz(-2.3402479) q[1];
sx q[1];
rz(-2.0739906) q[1];
sx q[1];
rz(1.9498922) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7758492) q[0];
sx q[0];
rz(-1.2677712) q[0];
sx q[0];
rz(-2.5822405) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0982047) q[2];
sx q[2];
rz(-2.049438) q[2];
sx q[2];
rz(-0.20028534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5687823) q[1];
sx q[1];
rz(-0.88925075) q[1];
sx q[1];
rz(3.0589025) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3016316) q[3];
sx q[3];
rz(-1.3084305) q[3];
sx q[3];
rz(0.056726168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2661065) q[2];
sx q[2];
rz(-2.2849639) q[2];
sx q[2];
rz(-0.53492707) q[2];
rz(2.5180425) q[3];
sx q[3];
rz(-2.0303191) q[3];
sx q[3];
rz(0.88304869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18735886) q[0];
sx q[0];
rz(-2.7281902) q[0];
sx q[0];
rz(0.9340539) q[0];
rz(-1.2454698) q[1];
sx q[1];
rz(-1.5555236) q[1];
sx q[1];
rz(-2.5159871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1067105) q[0];
sx q[0];
rz(-0.83256522) q[0];
sx q[0];
rz(-2.9009079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2410237) q[2];
sx q[2];
rz(-1.4975542) q[2];
sx q[2];
rz(1.1520958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.255068) q[1];
sx q[1];
rz(-2.2400948) q[1];
sx q[1];
rz(1.943735) q[1];
x q[2];
rz(0.91208338) q[3];
sx q[3];
rz(-2.8588223) q[3];
sx q[3];
rz(-1.8709489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0921359) q[2];
sx q[2];
rz(-2.226604) q[2];
sx q[2];
rz(-0.12506872) q[2];
rz(-0.72758979) q[3];
sx q[3];
rz(-1.5743419) q[3];
sx q[3];
rz(-2.6653813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8759988) q[0];
sx q[0];
rz(-2.6662874) q[0];
sx q[0];
rz(1.884961) q[0];
rz(-1.1988634) q[1];
sx q[1];
rz(-1.4224982) q[1];
sx q[1];
rz(-1.2748324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0755646) q[0];
sx q[0];
rz(-2.4504553) q[0];
sx q[0];
rz(0.71623556) q[0];
rz(-pi) q[1];
rz(2.149717) q[2];
sx q[2];
rz(-2.5235368) q[2];
sx q[2];
rz(1.2787873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9873571) q[1];
sx q[1];
rz(-1.2200095) q[1];
sx q[1];
rz(-1.1018176) q[1];
rz(1.2009462) q[3];
sx q[3];
rz(-2.8027034) q[3];
sx q[3];
rz(1.3047993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4510497) q[2];
sx q[2];
rz(-1.0131016) q[2];
sx q[2];
rz(0.96735442) q[2];
rz(1.5585772) q[3];
sx q[3];
rz(-1.6396921) q[3];
sx q[3];
rz(-0.30174747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3485182) q[0];
sx q[0];
rz(-0.59978849) q[0];
sx q[0];
rz(-0.10979688) q[0];
rz(1.9684567) q[1];
sx q[1];
rz(-2.1497612) q[1];
sx q[1];
rz(2.0593624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.968348) q[0];
sx q[0];
rz(-2.3254407) q[0];
sx q[0];
rz(1.9108755) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0556446) q[2];
sx q[2];
rz(-1.6669126) q[2];
sx q[2];
rz(-0.80730593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46690643) q[1];
sx q[1];
rz(-1.162858) q[1];
sx q[1];
rz(0.80927421) q[1];
x q[2];
rz(1.1940597) q[3];
sx q[3];
rz(-2.1749745) q[3];
sx q[3];
rz(0.24482152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3829696) q[2];
sx q[2];
rz(-1.0926282) q[2];
sx q[2];
rz(2.7782202) q[2];
rz(-2.162497) q[3];
sx q[3];
rz(-2.4978814) q[3];
sx q[3];
rz(2.6399829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30524224) q[0];
sx q[0];
rz(-2.3452121) q[0];
sx q[0];
rz(-2.7889732) q[0];
rz(1.7800219) q[1];
sx q[1];
rz(-1.1905328) q[1];
sx q[1];
rz(3.000066) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0240495) q[0];
sx q[0];
rz(-2.96857) q[0];
sx q[0];
rz(0.19917147) q[0];
rz(2.0317215) q[2];
sx q[2];
rz(-2.115098) q[2];
sx q[2];
rz(-1.5371295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6113811) q[1];
sx q[1];
rz(-0.57032835) q[1];
sx q[1];
rz(2.4228816) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64543313) q[3];
sx q[3];
rz(-1.107599) q[3];
sx q[3];
rz(2.2768159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8286459) q[2];
sx q[2];
rz(-1.6673648) q[2];
sx q[2];
rz(1.0860156) q[2];
rz(0.18276754) q[3];
sx q[3];
rz(-2.1153617) q[3];
sx q[3];
rz(-0.97203794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6163841) q[0];
sx q[0];
rz(-1.5807736) q[0];
sx q[0];
rz(-0.64424789) q[0];
rz(-3.0746025) q[1];
sx q[1];
rz(-1.3863775) q[1];
sx q[1];
rz(-0.11360528) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73175921) q[0];
sx q[0];
rz(-1.5749965) q[0];
sx q[0];
rz(-1.7948216) q[0];
rz(2.2502916) q[2];
sx q[2];
rz(-3.0105053) q[2];
sx q[2];
rz(-2.2830414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1190824) q[1];
sx q[1];
rz(-1.7412698) q[1];
sx q[1];
rz(-1.6891696) q[1];
x q[2];
rz(-1.9285517) q[3];
sx q[3];
rz(-1.610667) q[3];
sx q[3];
rz(1.8893582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99864787) q[2];
sx q[2];
rz(-1.8871658) q[2];
sx q[2];
rz(0.76510731) q[2];
rz(-2.9633925) q[3];
sx q[3];
rz(-1.5604115) q[3];
sx q[3];
rz(0.83711886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9255623) q[0];
sx q[0];
rz(-0.75556527) q[0];
sx q[0];
rz(2.8788993) q[0];
rz(-0.33949159) q[1];
sx q[1];
rz(-1.2295634) q[1];
sx q[1];
rz(-2.0801574) q[1];
rz(1.9148432) q[2];
sx q[2];
rz(-2.1969002) q[2];
sx q[2];
rz(0.67739396) q[2];
rz(-1.8858344) q[3];
sx q[3];
rz(-1.8585078) q[3];
sx q[3];
rz(0.97522492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
