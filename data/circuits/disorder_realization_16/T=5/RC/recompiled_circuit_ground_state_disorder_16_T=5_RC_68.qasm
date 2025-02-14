OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1333753) q[0];
sx q[0];
rz(-1.8741338) q[0];
sx q[0];
rz(-3.128669) q[0];
rz(-2.456993) q[1];
sx q[1];
rz(-0.79889387) q[1];
sx q[1];
rz(-1.0577143) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81257129) q[0];
sx q[0];
rz(-2.2264535) q[0];
sx q[0];
rz(-1.8581057) q[0];
x q[1];
rz(0.51527649) q[2];
sx q[2];
rz(-0.84723398) q[2];
sx q[2];
rz(-1.1267337) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2014241) q[1];
sx q[1];
rz(-0.88646171) q[1];
sx q[1];
rz(-2.0444872) q[1];
rz(1.6736044) q[3];
sx q[3];
rz(-0.3936201) q[3];
sx q[3];
rz(0.10904864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5554123) q[2];
sx q[2];
rz(-1.5988028) q[2];
sx q[2];
rz(-0.093322873) q[2];
rz(-3.1210476) q[3];
sx q[3];
rz(-2.8829657) q[3];
sx q[3];
rz(1.7790022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697295) q[0];
sx q[0];
rz(-1.7427895) q[0];
sx q[0];
rz(2.3139957) q[0];
rz(-0.96356511) q[1];
sx q[1];
rz(-1.6160485) q[1];
sx q[1];
rz(0.73659426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927836) q[0];
sx q[0];
rz(-0.28454706) q[0];
sx q[0];
rz(-1.297516) q[0];
x q[1];
rz(-0.97073769) q[2];
sx q[2];
rz(-2.063437) q[2];
sx q[2];
rz(-0.095183177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1087649) q[1];
sx q[1];
rz(-1.2695644) q[1];
sx q[1];
rz(2.9297622) q[1];
rz(-pi) q[2];
rz(-3.0071665) q[3];
sx q[3];
rz(-1.5930015) q[3];
sx q[3];
rz(-1.4060705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5654512) q[2];
sx q[2];
rz(-2.5665923) q[2];
sx q[2];
rz(-2.3973993) q[2];
rz(0.60892504) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(-1.5003381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610483) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(0.79958493) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(-2.0094357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47100859) q[0];
sx q[0];
rz(-0.23375227) q[0];
sx q[0];
rz(1.9712377) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8707809) q[2];
sx q[2];
rz(-0.083148227) q[2];
sx q[2];
rz(1.5025592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4375648) q[1];
sx q[1];
rz(-1.4933741) q[1];
sx q[1];
rz(-0.67005007) q[1];
rz(3.0789496) q[3];
sx q[3];
rz(-0.60986211) q[3];
sx q[3];
rz(2.8727639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75480294) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(1.7542138) q[2];
rz(-2.7925708) q[3];
sx q[3];
rz(-1.6882221) q[3];
sx q[3];
rz(-2.6121228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.741852) q[0];
sx q[0];
rz(-0.96804237) q[0];
sx q[0];
rz(-2.7834748) q[0];
rz(2.070836) q[1];
sx q[1];
rz(-2.4929969) q[1];
sx q[1];
rz(1.4395641) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69735295) q[0];
sx q[0];
rz(-1.3079738) q[0];
sx q[0];
rz(-0.79189827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1522572) q[2];
sx q[2];
rz(-2.3609567) q[2];
sx q[2];
rz(-0.5350185) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52709748) q[1];
sx q[1];
rz(-1.004809) q[1];
sx q[1];
rz(-2.4742592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.317512) q[3];
sx q[3];
rz(-2.321455) q[3];
sx q[3];
rz(1.2323605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4293356) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(-2.5168391) q[2];
rz(-1.5077) q[3];
sx q[3];
rz(-2.3816536) q[3];
sx q[3];
rz(-1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7655012) q[0];
sx q[0];
rz(-2.2608345) q[0];
sx q[0];
rz(-1.1464024) q[0];
rz(1.435185) q[1];
sx q[1];
rz(-2.621666) q[1];
sx q[1];
rz(-2.3211839) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56341237) q[0];
sx q[0];
rz(-0.2858361) q[0];
sx q[0];
rz(-1.8875627) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9632513) q[2];
sx q[2];
rz(-0.91195157) q[2];
sx q[2];
rz(2.4039011) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9238961) q[1];
sx q[1];
rz(-1.3753328) q[1];
sx q[1];
rz(0.20789897) q[1];
x q[2];
rz(0.049656258) q[3];
sx q[3];
rz(-0.59030246) q[3];
sx q[3];
rz(1.717339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4121805) q[2];
sx q[2];
rz(-1.0554487) q[2];
sx q[2];
rz(0.5385651) q[2];
rz(1.4696848) q[3];
sx q[3];
rz(-1.6939751) q[3];
sx q[3];
rz(0.058549747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3013714) q[0];
sx q[0];
rz(-0.138962) q[0];
sx q[0];
rz(-0.090601623) q[0];
rz(-0.36965707) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(-0.37818092) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8398741) q[0];
sx q[0];
rz(-1.0953971) q[0];
sx q[0];
rz(1.8283707) q[0];
x q[1];
rz(3.0514977) q[2];
sx q[2];
rz(-1.5290773) q[2];
sx q[2];
rz(0.99064529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82434618) q[1];
sx q[1];
rz(-0.50889824) q[1];
sx q[1];
rz(-2.8864278) q[1];
rz(-pi) q[2];
rz(3.0621594) q[3];
sx q[3];
rz(-1.7609247) q[3];
sx q[3];
rz(-3.0832689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52046627) q[2];
sx q[2];
rz(-1.4755321) q[2];
sx q[2];
rz(0.10144083) q[2];
rz(1.2927879) q[3];
sx q[3];
rz(-2.3088876) q[3];
sx q[3];
rz(-1.4147991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-1.3066278) q[0];
sx q[0];
rz(-1.8649626) q[0];
sx q[0];
rz(0.90245885) q[0];
rz(0.2392256) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(2.7801524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28736605) q[0];
sx q[0];
rz(-2.0660095) q[0];
sx q[0];
rz(-0.12676858) q[0];
x q[1];
rz(-3.0143033) q[2];
sx q[2];
rz(-1.5152928) q[2];
sx q[2];
rz(-1.0148359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4415008) q[1];
sx q[1];
rz(-0.48466408) q[1];
sx q[1];
rz(2.0551217) q[1];
rz(1.4124198) q[3];
sx q[3];
rz(-1.90622) q[3];
sx q[3];
rz(2.3349544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.07301894) q[2];
sx q[2];
rz(-1.8014427) q[2];
sx q[2];
rz(0.87693357) q[2];
rz(-0.98313156) q[3];
sx q[3];
rz(-1.5638331) q[3];
sx q[3];
rz(-2.1985998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.32901949) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(-2.3028497) q[0];
rz(0.70873952) q[1];
sx q[1];
rz(-2.510431) q[1];
sx q[1];
rz(-2.2355524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4137461) q[0];
sx q[0];
rz(-2.4165396) q[0];
sx q[0];
rz(-3.1097799) q[0];
rz(-pi) q[1];
rz(-2.110741) q[2];
sx q[2];
rz(-2.738852) q[2];
sx q[2];
rz(-2.7583721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4455394) q[1];
sx q[1];
rz(-1.0372682) q[1];
sx q[1];
rz(-0.75977709) q[1];
x q[2];
rz(0.25211199) q[3];
sx q[3];
rz(-2.1641755) q[3];
sx q[3];
rz(2.2098847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16443843) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(1.1158811) q[2];
rz(2.5130533) q[3];
sx q[3];
rz(-2.0597337) q[3];
sx q[3];
rz(2.5270497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81166613) q[0];
sx q[0];
rz(-0.8198494) q[0];
sx q[0];
rz(2.9803168) q[0];
rz(2.2992086) q[1];
sx q[1];
rz(-1.4295652) q[1];
sx q[1];
rz(-0.17002034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5710167) q[0];
sx q[0];
rz(-1.5525155) q[0];
sx q[0];
rz(1.5626426) q[0];
x q[1];
rz(2.7014637) q[2];
sx q[2];
rz(-1.7227731) q[2];
sx q[2];
rz(0.51896799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.57178088) q[1];
sx q[1];
rz(-1.7649635) q[1];
sx q[1];
rz(-2.012019) q[1];
rz(-1.1879425) q[3];
sx q[3];
rz(-1.597279) q[3];
sx q[3];
rz(-1.5838102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9876447) q[2];
sx q[2];
rz(-1.597007) q[2];
sx q[2];
rz(-1.6551931) q[2];
rz(-0.060128309) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.4415584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104367) q[0];
sx q[0];
rz(-3.1102409) q[0];
sx q[0];
rz(-1.7106868) q[0];
rz(-0.36262861) q[1];
sx q[1];
rz(-1.5321833) q[1];
sx q[1];
rz(1.3154202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.05689) q[0];
sx q[0];
rz(-1.2610208) q[0];
sx q[0];
rz(2.8357767) q[0];
x q[1];
rz(-1.270938) q[2];
sx q[2];
rz(-0.39405566) q[2];
sx q[2];
rz(1.9410417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23116701) q[1];
sx q[1];
rz(-2.1520414) q[1];
sx q[1];
rz(-1.2509996) q[1];
rz(0.4591367) q[3];
sx q[3];
rz(-0.92184421) q[3];
sx q[3];
rz(-2.5970363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.057283904) q[2];
sx q[2];
rz(-1.5949275) q[2];
sx q[2];
rz(2.9810737) q[2];
rz(-2.6594243) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(2.4619861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01874825) q[0];
sx q[0];
rz(-1.7457122) q[0];
sx q[0];
rz(2.2298298) q[0];
rz(1.1625166) q[1];
sx q[1];
rz(-1.569842) q[1];
sx q[1];
rz(-0.40649489) q[1];
rz(-0.087333655) q[2];
sx q[2];
rz(-1.4764016) q[2];
sx q[2];
rz(1.7511677) q[2];
rz(-1.0033458) q[3];
sx q[3];
rz(-0.88450817) q[3];
sx q[3];
rz(-1.3664288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
