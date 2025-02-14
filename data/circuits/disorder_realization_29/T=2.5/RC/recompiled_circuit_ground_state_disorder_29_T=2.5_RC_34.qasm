OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(5.7481264) q[0];
sx q[0];
rz(9.6587107) q[0];
rz(0.84199953) q[1];
sx q[1];
rz(-1.7276126) q[1];
sx q[1];
rz(1.6230621) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9568125) q[0];
sx q[0];
rz(-1.0597469) q[0];
sx q[0];
rz(1.2132116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0064924) q[2];
sx q[2];
rz(-2.6681136) q[2];
sx q[2];
rz(-1.4896637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7747353) q[1];
sx q[1];
rz(-2.1090057) q[1];
sx q[1];
rz(-2.4789206) q[1];
x q[2];
rz(-0.99125864) q[3];
sx q[3];
rz(-2.391444) q[3];
sx q[3];
rz(2.0969491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.4404094) q[2];
sx q[2];
rz(-2.3360628) q[2];
rz(-0.39189288) q[3];
sx q[3];
rz(-1.7854179) q[3];
sx q[3];
rz(-0.75683561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(-1.27502) q[1];
sx q[1];
rz(-1.6908815) q[1];
sx q[1];
rz(3.0335887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40168821) q[0];
sx q[0];
rz(-1.5459014) q[0];
sx q[0];
rz(-1.6817001) q[0];
rz(-pi) q[1];
rz(2.2453868) q[2];
sx q[2];
rz(-0.87088481) q[2];
sx q[2];
rz(2.7090379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9671772) q[1];
sx q[1];
rz(-1.6545873) q[1];
sx q[1];
rz(-1.3166974) q[1];
x q[2];
rz(-2.9505902) q[3];
sx q[3];
rz(-2.0536978) q[3];
sx q[3];
rz(-2.1306899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7763623) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(0.41415563) q[3];
sx q[3];
rz(-2.0678346) q[3];
sx q[3];
rz(-0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.47495833) q[0];
sx q[0];
rz(-1.8284429) q[0];
sx q[0];
rz(-2.3057002) q[0];
rz(-1.6235141) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(-1.2380884) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1956714) q[0];
sx q[0];
rz(-2.4950881) q[0];
sx q[0];
rz(-2.2158428) q[0];
rz(-pi) q[1];
rz(-2.5290501) q[2];
sx q[2];
rz(-0.41191891) q[2];
sx q[2];
rz(-1.4985794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1819396) q[1];
sx q[1];
rz(-1.6283855) q[1];
sx q[1];
rz(-1.4532386) q[1];
x q[2];
rz(-0.28210552) q[3];
sx q[3];
rz(-1.8477401) q[3];
sx q[3];
rz(-1.9953342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28430024) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(-0.71883744) q[2];
rz(-2.5607064) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(-0.7287997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.523664) q[0];
sx q[0];
rz(-1.5476462) q[0];
sx q[0];
rz(2.0148328) q[0];
rz(-1.2976546) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(2.0984971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422416) q[0];
sx q[0];
rz(-1.8411921) q[0];
sx q[0];
rz(-3.110136) q[0];
rz(-1.9773433) q[2];
sx q[2];
rz(-2.5669328) q[2];
sx q[2];
rz(2.3324147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14341893) q[1];
sx q[1];
rz(-0.34356782) q[1];
sx q[1];
rz(-2.0433389) q[1];
x q[2];
rz(-0.11123379) q[3];
sx q[3];
rz(-2.4832442) q[3];
sx q[3];
rz(1.2963795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2780693) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(-2.1045904) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(-0.83468848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0015513) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(-3.0539883) q[0];
rz(0.43830782) q[1];
sx q[1];
rz(-1.0821082) q[1];
sx q[1];
rz(1.8278488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3739819) q[0];
sx q[0];
rz(-1.0896519) q[0];
sx q[0];
rz(-0.26480459) q[0];
rz(1.472166) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(1.94869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.524595) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(1.3981546) q[1];
rz(-1.3145218) q[3];
sx q[3];
rz(-1.7328784) q[3];
sx q[3];
rz(2.7173998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4801415) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(0.56337774) q[2];
rz(1.2207458) q[3];
sx q[3];
rz(-1.3910339) q[3];
sx q[3];
rz(-2.5105072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0193943) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(-3.1296545) q[0];
rz(2.7519233) q[1];
sx q[1];
rz(-2.4395112) q[1];
sx q[1];
rz(0.24519244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58701506) q[0];
sx q[0];
rz(-2.3759807) q[0];
sx q[0];
rz(0.9258201) q[0];
x q[1];
rz(0.71227534) q[2];
sx q[2];
rz(-0.33604188) q[2];
sx q[2];
rz(-3.0672376) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17055146) q[1];
sx q[1];
rz(-2.1134796) q[1];
sx q[1];
rz(2.0849063) q[1];
rz(-2.7390476) q[3];
sx q[3];
rz(-1.8410826) q[3];
sx q[3];
rz(0.2337993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7616854) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(-1.3690108) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-2.4574418) q[3];
sx q[3];
rz(-2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928007) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(1.6311215) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(-2.7105601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1265701) q[0];
sx q[0];
rz(-1.3608785) q[0];
sx q[0];
rz(0.2951584) q[0];
rz(-1.6234342) q[2];
sx q[2];
rz(-1.8013012) q[2];
sx q[2];
rz(-0.33256691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37087155) q[1];
sx q[1];
rz(-1.8951192) q[1];
sx q[1];
rz(1.0350219) q[1];
x q[2];
rz(1.0930952) q[3];
sx q[3];
rz(-0.58148958) q[3];
sx q[3];
rz(0.026611004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83583528) q[2];
sx q[2];
rz(-0.77360669) q[2];
sx q[2];
rz(-0.26710278) q[2];
rz(0.032020656) q[3];
sx q[3];
rz(-1.9710385) q[3];
sx q[3];
rz(0.10572461) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06374643) q[0];
sx q[0];
rz(-1.8505322) q[0];
sx q[0];
rz(-2.5015976) q[0];
rz(-2.2867639) q[1];
sx q[1];
rz(-1.4850441) q[1];
sx q[1];
rz(1.7291501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873358) q[0];
sx q[0];
rz(-1.0414413) q[0];
sx q[0];
rz(-2.0482778) q[0];
rz(3.090254) q[2];
sx q[2];
rz(-1.079139) q[2];
sx q[2];
rz(1.3528388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.392726) q[1];
sx q[1];
rz(-0.35104529) q[1];
sx q[1];
rz(2.8944394) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7667609) q[3];
sx q[3];
rz(-1.9697947) q[3];
sx q[3];
rz(2.3712036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5414446) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(-2.7195462) q[2];
rz(-0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(0.1951018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37050978) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(-1.7899845) q[0];
rz(2.8260258) q[1];
sx q[1];
rz(-1.6866997) q[1];
sx q[1];
rz(1.9884761) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2446072) q[0];
sx q[0];
rz(-0.073052064) q[0];
sx q[0];
rz(0.79985072) q[0];
x q[1];
rz(-0.081772371) q[2];
sx q[2];
rz(-0.53682971) q[2];
sx q[2];
rz(-2.5935136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7866384) q[1];
sx q[1];
rz(-1.6980722) q[1];
sx q[1];
rz(2.5580225) q[1];
rz(-pi) q[2];
rz(-2.3117468) q[3];
sx q[3];
rz(-0.53433687) q[3];
sx q[3];
rz(-0.16610924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(0.37377629) q[2];
rz(-1.9155546) q[3];
sx q[3];
rz(-0.91807476) q[3];
sx q[3];
rz(1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1821197) q[0];
sx q[0];
rz(-2.4612893) q[0];
sx q[0];
rz(1.8208338) q[0];
rz(0.93718115) q[1];
sx q[1];
rz(-1.8056185) q[1];
sx q[1];
rz(-2.6722867) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.824911) q[0];
sx q[0];
rz(-1.4535722) q[0];
sx q[0];
rz(0.076923142) q[0];
rz(0.19005084) q[2];
sx q[2];
rz(-0.99744883) q[2];
sx q[2];
rz(-2.9721476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1782068) q[1];
sx q[1];
rz(-1.4904516) q[1];
sx q[1];
rz(-2.9909913) q[1];
rz(-0.68113459) q[3];
sx q[3];
rz(-2.7609842) q[3];
sx q[3];
rz(-1.8073624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4361973) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(-1.2104642) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.7567239) q[3];
sx q[3];
rz(-0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217011) q[0];
sx q[0];
rz(-0.85048631) q[0];
sx q[0];
rz(1.462107) q[0];
rz(2.2484491) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(-2.2371815) q[2];
sx q[2];
rz(-2.7795962) q[2];
sx q[2];
rz(-2.4301651) q[2];
rz(2.1410971) q[3];
sx q[3];
rz(-1.199493) q[3];
sx q[3];
rz(-1.5293157) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
