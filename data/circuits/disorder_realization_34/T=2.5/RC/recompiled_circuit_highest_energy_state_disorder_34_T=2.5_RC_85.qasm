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
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(0.42088977) q[0];
rz(-0.66840494) q[1];
sx q[1];
rz(4.6894046) q[1];
sx q[1];
rz(8.2735396) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1323038) q[0];
sx q[0];
rz(-1.124086) q[0];
sx q[0];
rz(-1.3473835) q[0];
x q[1];
rz(2.2368512) q[2];
sx q[2];
rz(-0.67585617) q[2];
sx q[2];
rz(-0.35206282) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80025339) q[1];
sx q[1];
rz(-1.0301939) q[1];
sx q[1];
rz(0.64250352) q[1];
x q[2];
rz(2.8731396) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(0.47273794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0737754) q[2];
sx q[2];
rz(-2.1442118) q[2];
sx q[2];
rz(-1.0214405) q[2];
rz(-1.8197618) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3385056) q[0];
sx q[0];
rz(-2.1144688) q[0];
sx q[0];
rz(-0.3013674) q[0];
rz(-2.001568) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(1.0637306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607644) q[0];
sx q[0];
rz(-1.6910292) q[0];
sx q[0];
rz(1.8390435) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3688886) q[2];
sx q[2];
rz(-0.81939745) q[2];
sx q[2];
rz(0.66172268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96460198) q[1];
sx q[1];
rz(-0.97787217) q[1];
sx q[1];
rz(-1.9511392) q[1];
rz(1.3724907) q[3];
sx q[3];
rz(-2.5970659) q[3];
sx q[3];
rz(2.209112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2782044) q[2];
sx q[2];
rz(-2.3953343) q[2];
sx q[2];
rz(-0.14872742) q[2];
rz(-1.9637828) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8389559) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(1.0021915) q[0];
rz(1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(-2.5340714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70456767) q[0];
sx q[0];
rz(-1.145073) q[0];
sx q[0];
rz(2.2680437) q[0];
rz(-pi) q[1];
rz(0.27670105) q[2];
sx q[2];
rz(-1.1778129) q[2];
sx q[2];
rz(0.0056841141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3866736) q[1];
sx q[1];
rz(-0.52820871) q[1];
sx q[1];
rz(-0.22294238) q[1];
rz(0.63849498) q[3];
sx q[3];
rz(-0.83052626) q[3];
sx q[3];
rz(2.3414224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9017631) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(0.98881161) q[2];
rz(1.6648939) q[3];
sx q[3];
rz(-1.8035382) q[3];
sx q[3];
rz(0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3636417) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(1.2203891) q[0];
rz(1.8266504) q[1];
sx q[1];
rz(-1.6362135) q[1];
sx q[1];
rz(3.0016518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264288) q[0];
sx q[0];
rz(-0.73472154) q[0];
sx q[0];
rz(0.61798851) q[0];
rz(-1.0331421) q[2];
sx q[2];
rz(-2.2058479) q[2];
sx q[2];
rz(0.5039353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9040457) q[1];
sx q[1];
rz(-1.1155459) q[1];
sx q[1];
rz(0.039467875) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071752944) q[3];
sx q[3];
rz(-1.9535616) q[3];
sx q[3];
rz(-2.0677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6949029) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(0.13738446) q[2];
rz(-1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3985652) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(-2.9803357) q[0];
rz(-0.84363168) q[1];
sx q[1];
rz(-2.3944941) q[1];
sx q[1];
rz(1.3074494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8358993) q[0];
sx q[0];
rz(-0.81863716) q[0];
sx q[0];
rz(0.37055117) q[0];
x q[1];
rz(2.9425386) q[2];
sx q[2];
rz(-1.2110333) q[2];
sx q[2];
rz(-2.2729286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89462639) q[1];
sx q[1];
rz(-1.6603396) q[1];
sx q[1];
rz(0.98693165) q[1];
x q[2];
rz(0.90353538) q[3];
sx q[3];
rz(-1.0097111) q[3];
sx q[3];
rz(-0.72431952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1696986) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(0.92740721) q[2];
rz(1.6124604) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(-2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2138432) q[0];
sx q[0];
rz(-0.96927154) q[0];
sx q[0];
rz(1.4056322) q[0];
rz(-1.4145781) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(-2.3419211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7223631) q[0];
sx q[0];
rz(-1.5088313) q[0];
sx q[0];
rz(-0.011716066) q[0];
rz(-0.39787103) q[2];
sx q[2];
rz(-2.8074773) q[2];
sx q[2];
rz(-1.0256922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1649014) q[1];
sx q[1];
rz(-1.8215239) q[1];
sx q[1];
rz(-1.2198971) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7148083) q[3];
sx q[3];
rz(-0.43922654) q[3];
sx q[3];
rz(2.3595485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8863525) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(0.8026455) q[2];
rz(-2.8211527) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(-1.6360487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36444148) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(0.030315422) q[0];
rz(1.8346571) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(2.511715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4277074) q[0];
sx q[0];
rz(-0.46525912) q[0];
sx q[0];
rz(2.7119539) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3223619) q[2];
sx q[2];
rz(-2.5360245) q[2];
sx q[2];
rz(1.8964963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6678289) q[1];
sx q[1];
rz(-2.7817717) q[1];
sx q[1];
rz(-2.5303165) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23655741) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(2.5270324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1097172) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(-1.1472222) q[2];
rz(-2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77618337) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(2.6614406) q[0];
rz(-2.7393553) q[1];
sx q[1];
rz(-1.4996585) q[1];
sx q[1];
rz(2.7587836) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881809) q[0];
sx q[0];
rz(-2.7482902) q[0];
sx q[0];
rz(-2.1476782) q[0];
x q[1];
rz(0.11522861) q[2];
sx q[2];
rz(-1.281321) q[2];
sx q[2];
rz(2.8897499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9342707) q[1];
sx q[1];
rz(-1.2756961) q[1];
sx q[1];
rz(1.9701824) q[1];
rz(1.6005101) q[3];
sx q[3];
rz(-1.7669356) q[3];
sx q[3];
rz(0.33523424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.39168921) q[2];
sx q[2];
rz(-2.4112371) q[2];
sx q[2];
rz(-0.35823092) q[2];
rz(0.41424888) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888032) q[0];
sx q[0];
rz(-1.8599334) q[0];
sx q[0];
rz(-2.7237256) q[0];
rz(-1.4990384) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(-0.61161673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30464028) q[0];
sx q[0];
rz(-1.3654183) q[0];
sx q[0];
rz(-2.9083328) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29199227) q[2];
sx q[2];
rz(-2.9270009) q[2];
sx q[2];
rz(0.24927441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4222718) q[1];
sx q[1];
rz(-0.55972717) q[1];
sx q[1];
rz(-2.1156963) q[1];
rz(-2.9183594) q[3];
sx q[3];
rz(-2.7866552) q[3];
sx q[3];
rz(2.8508773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.186782) q[2];
sx q[2];
rz(-1.45767) q[2];
sx q[2];
rz(0.35624722) q[2];
rz(-0.48480836) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(-0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.818882) q[0];
sx q[0];
rz(-2.0368545) q[0];
sx q[0];
rz(-1.6792962) q[0];
rz(-2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(-0.80518728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96435409) q[0];
sx q[0];
rz(-0.88429175) q[0];
sx q[0];
rz(0.80857386) q[0];
rz(-0.83909713) q[2];
sx q[2];
rz(-0.57575127) q[2];
sx q[2];
rz(-2.8481399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78811121) q[1];
sx q[1];
rz(-0.32057163) q[1];
sx q[1];
rz(0.48133565) q[1];
x q[2];
rz(-2.7789742) q[3];
sx q[3];
rz(-1.3380074) q[3];
sx q[3];
rz(-1.5513368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8055973) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(-1.9192637) q[2];
rz(-0.81021106) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(1.5626102) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41351086) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(-1.8660846) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(-0.13678094) q[2];
sx q[2];
rz(-1.9521458) q[2];
sx q[2];
rz(-1.4804179) q[2];
rz(-2.2386407) q[3];
sx q[3];
rz(-1.9330238) q[3];
sx q[3];
rz(-0.15180363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
