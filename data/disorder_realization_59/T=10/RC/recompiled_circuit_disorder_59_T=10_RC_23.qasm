OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(0.0013874887) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35933094) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(0.1749392) q[0];
rz(-pi) q[1];
rz(-0.44714655) q[2];
sx q[2];
rz(-1.5329754) q[2];
sx q[2];
rz(0.58085261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5056155) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(-2.99519) q[1];
rz(0.9790768) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(0.81545365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0464697) q[0];
sx q[0];
rz(-0.9233343) q[0];
sx q[0];
rz(1.6460653) q[0];
x q[1];
rz(-0.30351992) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(2.7941861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26317877) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(3.0221699) q[1];
x q[2];
rz(0.64579441) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(-1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(-1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777269) q[0];
sx q[0];
rz(-2.8054872) q[0];
sx q[0];
rz(1.1019812) q[0];
rz(-0.83061647) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-2.9142771) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(-1.2815777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0619034) q[3];
sx q[3];
rz(-1.5012) q[3];
sx q[3];
rz(0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-0.27483637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1823605) q[0];
sx q[0];
rz(-1.649547) q[0];
sx q[0];
rz(1.4819281) q[0];
rz(-pi) q[1];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(-1.9542076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6064925) q[1];
sx q[1];
rz(-0.98787687) q[1];
sx q[1];
rz(0.54775723) q[1];
x q[2];
rz(-1.4196017) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.5083195) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65790025) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18407962) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(-2.1225131) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6078569) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.5460207) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21302528) q[1];
sx q[1];
rz(-1.2270317) q[1];
sx q[1];
rz(-2.7108971) q[1];
rz(0.98832163) q[3];
sx q[3];
rz(-2.0867996) q[3];
sx q[3];
rz(0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-2.8175763) q[2];
rz(1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(-0.34067672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5013803) q[0];
sx q[0];
rz(-0.073177241) q[0];
sx q[0];
rz(2.0206547) q[0];
x q[1];
rz(-0.31008115) q[2];
sx q[2];
rz(-0.78044621) q[2];
sx q[2];
rz(0.23852894) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1250455) q[1];
sx q[1];
rz(-2.6066337) q[1];
sx q[1];
rz(-0.46334456) q[1];
rz(-2.0601222) q[3];
sx q[3];
rz(-1.9483856) q[3];
sx q[3];
rz(1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(-2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-0.46494928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5349605) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(1.1688933) q[0];
rz(-pi) q[1];
rz(2.6480688) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.354419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.773015) q[1];
sx q[1];
rz(-1.7891208) q[1];
sx q[1];
rz(1.6792084) q[1];
rz(-0.025860272) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6440455) q[0];
sx q[0];
rz(-0.80701485) q[0];
sx q[0];
rz(3.0456196) q[0];
x q[1];
rz(2.3381091) q[2];
sx q[2];
rz(-2.8065971) q[2];
sx q[2];
rz(1.6840881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6137177) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(-1.2916958) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(0.3233288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.4096227) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(0.87402469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4955935) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(1.5587224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2888695) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(-0.66327099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28305001) q[1];
sx q[1];
rz(-1.0714604) q[1];
sx q[1];
rz(0.081375558) q[1];
rz(-pi) q[2];
rz(2.5328818) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.727227) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.1788517) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26319474) q[0];
sx q[0];
rz(-0.34391719) q[0];
sx q[0];
rz(-0.92349903) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6640501) q[2];
sx q[2];
rz(-1.6972731) q[2];
sx q[2];
rz(1.3716413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5621592) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.8821554) q[1];
rz(1.3445271) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(-2.1567878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205004) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.6451251) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
