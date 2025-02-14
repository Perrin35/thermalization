OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3235029) q[0];
sx q[0];
rz(-2.7348195) q[0];
sx q[0];
rz(2.6320501) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(-2.9480313) q[1];
sx q[1];
rz(-2.136018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5339342) q[0];
sx q[0];
rz(-2.9276121) q[0];
sx q[0];
rz(-0.35983742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2174913) q[2];
sx q[2];
rz(-1.2900724) q[2];
sx q[2];
rz(-1.4626056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7889346) q[1];
sx q[1];
rz(-0.96436687) q[1];
sx q[1];
rz(-2.2878617) q[1];
x q[2];
rz(-2.5343603) q[3];
sx q[3];
rz(-2.2650026) q[3];
sx q[3];
rz(1.1986365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84844184) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(-1.1085917) q[2];
rz(-2.0075924) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(3.0602684) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14652458) q[0];
sx q[0];
rz(-1.1256555) q[0];
sx q[0];
rz(2.8296237) q[0];
rz(0.23090714) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(-1.1280967) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8618516) q[0];
sx q[0];
rz(-1.3817046) q[0];
sx q[0];
rz(-1.0936589) q[0];
rz(0.53324576) q[2];
sx q[2];
rz(-0.74138481) q[2];
sx q[2];
rz(-0.35517755) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1027185) q[1];
sx q[1];
rz(-0.46157122) q[1];
sx q[1];
rz(-0.25074236) q[1];
x q[2];
rz(0.76448582) q[3];
sx q[3];
rz(-1.8967046) q[3];
sx q[3];
rz(-0.21316646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6850623) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(-0.064229639) q[2];
rz(1.8188933) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(-0.88671154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.54881683) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(-1.892426) q[0];
rz(0.33117548) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(-1.9452728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0421233) q[0];
sx q[0];
rz(-1.5705918) q[0];
sx q[0];
rz(1.5705137) q[0];
x q[1];
rz(1.565479) q[2];
sx q[2];
rz(-1.8182767) q[2];
sx q[2];
rz(0.70308816) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0532562) q[1];
sx q[1];
rz(-2.1532946) q[1];
sx q[1];
rz(0.94515578) q[1];
rz(2.7306741) q[3];
sx q[3];
rz(-2.344729) q[3];
sx q[3];
rz(-1.381402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(-2.2288442) q[2];
rz(1.4384455) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(1.335817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790633) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(-0.67034876) q[0];
rz(2.4743075) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(0.60428062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33907783) q[0];
sx q[0];
rz(-2.3065563) q[0];
sx q[0];
rz(-1.5123532) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3717643) q[2];
sx q[2];
rz(-0.40707874) q[2];
sx q[2];
rz(2.4224015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5126219) q[1];
sx q[1];
rz(-1.781946) q[1];
sx q[1];
rz(1.8308012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8205582) q[3];
sx q[3];
rz(-0.74640025) q[3];
sx q[3];
rz(1.6370893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1112572) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(0.37720171) q[2];
rz(-3.0180569) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443026) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(-2.8422728) q[0];
rz(-2.0206644) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(2.1790806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6448005) q[0];
sx q[0];
rz(-1.3649277) q[0];
sx q[0];
rz(1.0394761) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6766432) q[2];
sx q[2];
rz(-2.1735149) q[2];
sx q[2];
rz(-0.51747396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.968041) q[1];
sx q[1];
rz(-2.3690201) q[1];
sx q[1];
rz(-2.3501758) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59904091) q[3];
sx q[3];
rz(-1.1223464) q[3];
sx q[3];
rz(0.34235172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0472497) q[2];
sx q[2];
rz(-0.80545682) q[2];
sx q[2];
rz(-0.94364014) q[2];
rz(1.0654248) q[3];
sx q[3];
rz(-1.7908955) q[3];
sx q[3];
rz(-2.0431199) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0672673) q[0];
sx q[0];
rz(-1.4370947) q[0];
sx q[0];
rz(3.1332916) q[0];
rz(0.51757327) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(1.5483206) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8529274) q[0];
sx q[0];
rz(-1.2199243) q[0];
sx q[0];
rz(-0.025625833) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6363434) q[2];
sx q[2];
rz(-1.3224241) q[2];
sx q[2];
rz(-0.58282436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0071809) q[1];
sx q[1];
rz(-2.6406277) q[1];
sx q[1];
rz(-2.3208614) q[1];
rz(-pi) q[2];
rz(0.41554873) q[3];
sx q[3];
rz(-2.1449148) q[3];
sx q[3];
rz(-2.1972063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3090618) q[2];
sx q[2];
rz(-1.4893724) q[2];
sx q[2];
rz(1.8050516) q[2];
rz(3.0858223) q[3];
sx q[3];
rz(-2.1499108) q[3];
sx q[3];
rz(-0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5830773) q[0];
sx q[0];
rz(-2.6641088) q[0];
sx q[0];
rz(0.68786311) q[0];
rz(-1.6784003) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(-2.8942143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32175049) q[0];
sx q[0];
rz(-0.30748707) q[0];
sx q[0];
rz(-2.5006258) q[0];
rz(-pi) q[1];
x q[1];
rz(2.192657) q[2];
sx q[2];
rz(-1.0100216) q[2];
sx q[2];
rz(-1.6552629) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.069889594) q[1];
sx q[1];
rz(-0.99537288) q[1];
sx q[1];
rz(-2.6399122) q[1];
x q[2];
rz(-1.3379657) q[3];
sx q[3];
rz(-1.4484753) q[3];
sx q[3];
rz(1.2965152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6972202) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(-0.42416254) q[2];
rz(-1.0423202) q[3];
sx q[3];
rz(-2.3764231) q[3];
sx q[3];
rz(-0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384077) q[0];
sx q[0];
rz(-0.98588949) q[0];
sx q[0];
rz(0.61088046) q[0];
rz(-0.25018397) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(0.47666034) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7358897) q[0];
sx q[0];
rz(-0.48375706) q[0];
sx q[0];
rz(1.4866845) q[0];
rz(-pi) q[1];
rz(-1.119461) q[2];
sx q[2];
rz(-2.9381466) q[2];
sx q[2];
rz(2.1687393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3150419) q[1];
sx q[1];
rz(-1.0796865) q[1];
sx q[1];
rz(1.1129614) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42929452) q[3];
sx q[3];
rz(-1.0280812) q[3];
sx q[3];
rz(2.8739704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99159795) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-2.4617713) q[2];
rz(-0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(1.6405039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221472) q[0];
sx q[0];
rz(-1.7663904) q[0];
sx q[0];
rz(2.9875901) q[0];
rz(0.18140659) q[1];
sx q[1];
rz(-1.0702952) q[1];
sx q[1];
rz(3.0214686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4331872) q[0];
sx q[0];
rz(-0.4077724) q[0];
sx q[0];
rz(-1.4207409) q[0];
rz(1.3072877) q[2];
sx q[2];
rz(-1.3804091) q[2];
sx q[2];
rz(0.15979494) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57615818) q[1];
sx q[1];
rz(-1.6685969) q[1];
sx q[1];
rz(1.3937772) q[1];
x q[2];
rz(-1.7150015) q[3];
sx q[3];
rz(-2.0940082) q[3];
sx q[3];
rz(-1.6169523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7141562) q[2];
sx q[2];
rz(-1.7566046) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(-0.51042026) q[3];
sx q[3];
rz(-2.3000058) q[3];
sx q[3];
rz(0.69971219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0521312) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(-2.7556038) q[0];
rz(1.5486859) q[1];
sx q[1];
rz(-0.49647757) q[1];
sx q[1];
rz(1.5923502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95741877) q[0];
sx q[0];
rz(-0.17137073) q[0];
sx q[0];
rz(1.4032073) q[0];
x q[1];
rz(-0.73076325) q[2];
sx q[2];
rz(-0.77789069) q[2];
sx q[2];
rz(0.90532263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4262271) q[1];
sx q[1];
rz(-1.5439529) q[1];
sx q[1];
rz(-1.6728841) q[1];
rz(-1.9796124) q[3];
sx q[3];
rz(-2.9745515) q[3];
sx q[3];
rz(0.80834889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3647032) q[2];
sx q[2];
rz(-0.93901712) q[2];
sx q[2];
rz(-1.6949863) q[2];
rz(1.0363091) q[3];
sx q[3];
rz(-2.2615137) q[3];
sx q[3];
rz(-0.026329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25630367) q[0];
sx q[0];
rz(-2.1623609) q[0];
sx q[0];
rz(1.4872861) q[0];
rz(-2.9215095) q[1];
sx q[1];
rz(-2.0368123) q[1];
sx q[1];
rz(-1.4727551) q[1];
rz(0.73843602) q[2];
sx q[2];
rz(-1.8059352) q[2];
sx q[2];
rz(-0.21485938) q[2];
rz(-1.2418048) q[3];
sx q[3];
rz(-2.0488779) q[3];
sx q[3];
rz(-2.4490801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
