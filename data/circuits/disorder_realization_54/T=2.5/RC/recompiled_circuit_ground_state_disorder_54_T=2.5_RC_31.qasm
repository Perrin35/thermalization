OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(5.1533617) q[0];
sx q[0];
rz(12.552153) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(0.74906936) q[1];
sx q[1];
rz(11.963117) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8612807) q[0];
sx q[0];
rz(-0.93398184) q[0];
sx q[0];
rz(2.0122819) q[0];
rz(-pi) q[1];
rz(1.649734) q[2];
sx q[2];
rz(-2.4853155) q[2];
sx q[2];
rz(0.62103926) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95663623) q[1];
sx q[1];
rz(-2.2549964) q[1];
sx q[1];
rz(1.7723899) q[1];
rz(-pi) q[2];
rz(-1.7885277) q[3];
sx q[3];
rz(-0.7184295) q[3];
sx q[3];
rz(-1.4361385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3732036) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(1.0962037) q[2];
rz(-2.0632035) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(-0.53511846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815149) q[0];
sx q[0];
rz(-2.09477) q[0];
sx q[0];
rz(2.6640418) q[0];
rz(2.1304456) q[1];
sx q[1];
rz(-0.54713455) q[1];
sx q[1];
rz(-2.0571041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80637425) q[0];
sx q[0];
rz(-1.9382297) q[0];
sx q[0];
rz(3.0159877) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5039551) q[2];
sx q[2];
rz(-2.7324841) q[2];
sx q[2];
rz(0.44873777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66531679) q[1];
sx q[1];
rz(-1.3937104) q[1];
sx q[1];
rz(2.5873273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73954496) q[3];
sx q[3];
rz(-1.7266183) q[3];
sx q[3];
rz(-0.23529822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61887211) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(1.0901394) q[2];
rz(-2.5777396) q[3];
sx q[3];
rz(-2.4520051) q[3];
sx q[3];
rz(-3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801333) q[0];
sx q[0];
rz(-3.1205803) q[0];
sx q[0];
rz(-1.6047961) q[0];
rz(1.3179294) q[1];
sx q[1];
rz(-1.0937966) q[1];
sx q[1];
rz(-0.39753786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73765342) q[0];
sx q[0];
rz(-1.353147) q[0];
sx q[0];
rz(-2.2277839) q[0];
x q[1];
rz(-0.99293844) q[2];
sx q[2];
rz(-1.8210337) q[2];
sx q[2];
rz(1.1838352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2469625) q[1];
sx q[1];
rz(-1.3161462) q[1];
sx q[1];
rz(0.92856426) q[1];
x q[2];
rz(-0.52513772) q[3];
sx q[3];
rz(-2.0179837) q[3];
sx q[3];
rz(-0.95765985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90618187) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(-2.4681674) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(1.5448236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509318) q[0];
sx q[0];
rz(-1.1216102) q[0];
sx q[0];
rz(2.2326873) q[0];
rz(0.017223651) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-3.1294894) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41866583) q[0];
sx q[0];
rz(-1.2483435) q[0];
sx q[0];
rz(-2.5315383) q[0];
rz(-pi) q[1];
rz(-0.64938371) q[2];
sx q[2];
rz(-2.0133791) q[2];
sx q[2];
rz(-1.8049406) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3693272) q[1];
sx q[1];
rz(-0.64913926) q[1];
sx q[1];
rz(-1.1916129) q[1];
rz(3.0148325) q[3];
sx q[3];
rz(-1.6157627) q[3];
sx q[3];
rz(0.25925207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072731344) q[2];
sx q[2];
rz(-2.8673745) q[2];
sx q[2];
rz(0.38632986) q[2];
rz(2.7522411) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(-2.1336335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35578457) q[0];
sx q[0];
rz(-0.72440994) q[0];
sx q[0];
rz(2.8177148) q[0];
rz(2.7549699) q[1];
sx q[1];
rz(-1.3893678) q[1];
sx q[1];
rz(1.5547543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51832685) q[0];
sx q[0];
rz(-2.4072106) q[0];
sx q[0];
rz(0.096303864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0350758) q[2];
sx q[2];
rz(-1.8958099) q[2];
sx q[2];
rz(-2.6485788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.388098) q[1];
sx q[1];
rz(-2.3109204) q[1];
sx q[1];
rz(2.4945033) q[1];
rz(-pi) q[2];
rz(-0.50916785) q[3];
sx q[3];
rz(-2.2751791) q[3];
sx q[3];
rz(-0.99607638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(-3.0774934) q[2];
rz(0.60424232) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(-1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20246777) q[0];
sx q[0];
rz(-0.61102837) q[0];
sx q[0];
rz(2.2623999) q[0];
rz(2.7912256) q[1];
sx q[1];
rz(-1.3060952) q[1];
sx q[1];
rz(2.9894357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504753) q[0];
sx q[0];
rz(-0.88248173) q[0];
sx q[0];
rz(1.2173714) q[0];
rz(0.19292508) q[2];
sx q[2];
rz(-1.0954276) q[2];
sx q[2];
rz(2.1383291) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2895246) q[1];
sx q[1];
rz(-1.7729323) q[1];
sx q[1];
rz(-0.44092559) q[1];
x q[2];
rz(-0.44701266) q[3];
sx q[3];
rz(-1.2114297) q[3];
sx q[3];
rz(-1.4554909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0396314) q[2];
sx q[2];
rz(-0.32932082) q[2];
sx q[2];
rz(2.9015818) q[2];
rz(2.4890684) q[3];
sx q[3];
rz(-1.1381166) q[3];
sx q[3];
rz(1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(-2.2834593) q[0];
rz(-1.7543322) q[1];
sx q[1];
rz(-2.6871082) q[1];
sx q[1];
rz(0.33285704) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8008566) q[0];
sx q[0];
rz(-1.256167) q[0];
sx q[0];
rz(2.1701909) q[0];
x q[1];
rz(-0.22631876) q[2];
sx q[2];
rz(-2.8853325) q[2];
sx q[2];
rz(-1.2761436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.27204) q[1];
sx q[1];
rz(-1.7764047) q[1];
sx q[1];
rz(2.7133792) q[1];
rz(-pi) q[2];
rz(2.694533) q[3];
sx q[3];
rz(-0.70526988) q[3];
sx q[3];
rz(1.3652319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.7994999) q[2];
sx q[2];
rz(-0.19793333) q[2];
sx q[2];
rz(-2.0905154) q[2];
rz(-1.6445232) q[3];
sx q[3];
rz(-1.1727419) q[3];
sx q[3];
rz(0.79819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.368211) q[0];
sx q[0];
rz(-1.0268651) q[0];
sx q[0];
rz(-2.2494702) q[0];
rz(0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(0.77286744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65281212) q[0];
sx q[0];
rz(-2.1202181) q[0];
sx q[0];
rz(2.1958952) q[0];
rz(1.8310089) q[2];
sx q[2];
rz(-0.88426829) q[2];
sx q[2];
rz(1.7169184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1972496) q[1];
sx q[1];
rz(-2.5641003) q[1];
sx q[1];
rz(0.61402278) q[1];
rz(-2.3691133) q[3];
sx q[3];
rz(-2.8101343) q[3];
sx q[3];
rz(0.14642388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.2246776) q[2];
sx q[2];
rz(-0.69250715) q[2];
rz(-0.45037371) q[3];
sx q[3];
rz(-1.2053442) q[3];
sx q[3];
rz(-1.5041806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5699128) q[0];
sx q[0];
rz(-0.79140651) q[0];
sx q[0];
rz(0.94863844) q[0];
rz(-1.0665077) q[1];
sx q[1];
rz(-0.98943168) q[1];
sx q[1];
rz(1.271064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9635634) q[0];
sx q[0];
rz(-0.22582136) q[0];
sx q[0];
rz(-1.9972576) q[0];
rz(1.5367212) q[2];
sx q[2];
rz(-1.5487031) q[2];
sx q[2];
rz(1.3349229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.003351) q[1];
sx q[1];
rz(-1.4356704) q[1];
sx q[1];
rz(-2.1993932) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27008121) q[3];
sx q[3];
rz(-2.684786) q[3];
sx q[3];
rz(-2.1169259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5991685) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(3.001281) q[2];
rz(3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(-1.9353346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5133544) q[0];
sx q[0];
rz(-1.6980549) q[0];
sx q[0];
rz(0.7775318) q[0];
rz(-0.93742049) q[1];
sx q[1];
rz(-1.9098858) q[1];
sx q[1];
rz(-2.6984528) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4454058) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(0.085025351) q[0];
rz(-pi) q[1];
rz(0.97340092) q[2];
sx q[2];
rz(-1.294916) q[2];
sx q[2];
rz(1.2593002) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75966802) q[1];
sx q[1];
rz(-0.72205359) q[1];
sx q[1];
rz(-0.51314919) q[1];
rz(-pi) q[2];
rz(2.0432161) q[3];
sx q[3];
rz(-2.7948423) q[3];
sx q[3];
rz(0.43671331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-3.013986) q[2];
rz(2.8429032) q[3];
sx q[3];
rz(-1.1726215) q[3];
sx q[3];
rz(0.3704139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5304607) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(-2.5313189) q[1];
sx q[1];
rz(-2.2563969) q[1];
sx q[1];
rz(1.0856249) q[1];
rz(2.2895428) q[2];
sx q[2];
rz(-1.5222737) q[2];
sx q[2];
rz(0.52494502) q[2];
rz(0.21214938) q[3];
sx q[3];
rz(-1.129088) q[3];
sx q[3];
rz(-1.3442985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
