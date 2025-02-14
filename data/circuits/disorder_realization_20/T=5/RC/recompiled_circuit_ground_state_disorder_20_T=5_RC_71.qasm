OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(0.33302745) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(3.9685213) q[1];
sx q[1];
rz(10.068738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4304023) q[0];
sx q[0];
rz(-1.2997753) q[0];
sx q[0];
rz(-0.94078101) q[0];
rz(-2.8415235) q[2];
sx q[2];
rz(-1.7207533) q[2];
sx q[2];
rz(-1.7008925) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17929303) q[1];
sx q[1];
rz(-2.6897618) q[1];
sx q[1];
rz(-2.5746114) q[1];
rz(-0.87857492) q[3];
sx q[3];
rz(-0.51673104) q[3];
sx q[3];
rz(-1.0754835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6551299) q[2];
sx q[2];
rz(-0.4643521) q[2];
sx q[2];
rz(-0.74074024) q[2];
rz(1.5898534) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(-2.5488502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43316677) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(-3.0246227) q[0];
rz(0.50432694) q[1];
sx q[1];
rz(-1.3386936) q[1];
sx q[1];
rz(0.28409827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0237673) q[0];
sx q[0];
rz(-1.6540915) q[0];
sx q[0];
rz(2.1375594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1524736) q[2];
sx q[2];
rz(-2.2788413) q[2];
sx q[2];
rz(3.0281554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71775466) q[1];
sx q[1];
rz(-1.6053797) q[1];
sx q[1];
rz(1.6929021) q[1];
x q[2];
rz(0.25155622) q[3];
sx q[3];
rz(-1.8402035) q[3];
sx q[3];
rz(3.0602853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8154907) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(-0.76078129) q[2];
rz(2.271999) q[3];
sx q[3];
rz(-1.2660916) q[3];
sx q[3];
rz(-2.6884955) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3608383) q[0];
sx q[0];
rz(-1.1203082) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(1.4422656) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(-0.35626492) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748791) q[0];
sx q[0];
rz(-1.5646114) q[0];
sx q[0];
rz(-1.6303819) q[0];
rz(2.5222657) q[2];
sx q[2];
rz(-2.4973923) q[2];
sx q[2];
rz(-1.2540115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9911451) q[1];
sx q[1];
rz(-0.65288645) q[1];
sx q[1];
rz(-1.3351403) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3144794) q[3];
sx q[3];
rz(-1.5689092) q[3];
sx q[3];
rz(-2.6607115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0226655) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(-1.5125795) q[2];
rz(3.1318943) q[3];
sx q[3];
rz(-1.2303979) q[3];
sx q[3];
rz(2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951204) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(2.8908308) q[0];
rz(1.7950902) q[1];
sx q[1];
rz(-2.612412) q[1];
sx q[1];
rz(1.4432602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478717) q[0];
sx q[0];
rz(-1.3547055) q[0];
sx q[0];
rz(-2.3225075) q[0];
rz(0.40374229) q[2];
sx q[2];
rz(-1.0761257) q[2];
sx q[2];
rz(0.97569377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6175198) q[1];
sx q[1];
rz(-2.1012602) q[1];
sx q[1];
rz(0.73442119) q[1];
x q[2];
rz(-2.3425927) q[3];
sx q[3];
rz(-2.7741787) q[3];
sx q[3];
rz(-2.011881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5556339) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(-2.8982758) q[2];
rz(0.58468741) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(0.028133597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(2.962033) q[0];
rz(-1.2252294) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(2.6833351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97044635) q[0];
sx q[0];
rz(-1.866328) q[0];
sx q[0];
rz(-1.4895205) q[0];
x q[1];
rz(2.8446537) q[2];
sx q[2];
rz(-0.7504979) q[2];
sx q[2];
rz(0.82210449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0441664) q[1];
sx q[1];
rz(-1.4961745) q[1];
sx q[1];
rz(1.7057555) q[1];
rz(-pi) q[2];
rz(2.7814976) q[3];
sx q[3];
rz(-2.2392139) q[3];
sx q[3];
rz(2.6289712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7563584) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(0.9217841) q[2];
rz(1.8900169) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0475622) q[0];
sx q[0];
rz(-2.229409) q[0];
sx q[0];
rz(-2.9929274) q[0];
rz(-0.31691638) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(-2.5247578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6233404) q[0];
sx q[0];
rz(-0.099661544) q[0];
sx q[0];
rz(-0.90604337) q[0];
rz(1.1304752) q[2];
sx q[2];
rz(-1.2424349) q[2];
sx q[2];
rz(-2.2715914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1779536) q[1];
sx q[1];
rz(-2.511189) q[1];
sx q[1];
rz(1.1373345) q[1];
x q[2];
rz(-2.0265977) q[3];
sx q[3];
rz(-0.57145703) q[3];
sx q[3];
rz(0.28459099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(3.1332664) q[2];
rz(2.5152123) q[3];
sx q[3];
rz(-0.71599394) q[3];
sx q[3];
rz(0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.9897113) q[0];
sx q[0];
rz(-2.6595111) q[0];
rz(-3.112402) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(0.76404244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.5034165) q[0];
sx q[0];
rz(-1.9420366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2305626) q[2];
sx q[2];
rz(-1.5907968) q[2];
sx q[2];
rz(2.3338712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58043874) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.2777722) q[1];
rz(-3.1336083) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(-0.29618357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6782516) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(-0.48416644) q[2];
rz(0.68228996) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(0.13474034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(-2.6765587) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(2.8344287) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7028593) q[0];
sx q[0];
rz(-0.68638681) q[0];
sx q[0];
rz(-2.2728361) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1910153) q[2];
sx q[2];
rz(-0.75846106) q[2];
sx q[2];
rz(-1.72067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0867566) q[1];
sx q[1];
rz(-2.1385178) q[1];
sx q[1];
rz(0.18594976) q[1];
x q[2];
rz(0.028178111) q[3];
sx q[3];
rz(-2.5163076) q[3];
sx q[3];
rz(0.5042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18053599) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(0.72009909) q[2];
rz(0.85047203) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(1.7299962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9367323) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(1.5180961) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(-2.785397) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30763876) q[0];
sx q[0];
rz(-2.4705268) q[0];
sx q[0];
rz(2.4249581) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61788606) q[2];
sx q[2];
rz(-1.5949342) q[2];
sx q[2];
rz(-1.2302421) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56936383) q[1];
sx q[1];
rz(-2.4920643) q[1];
sx q[1];
rz(0.049687548) q[1];
rz(1.1061927) q[3];
sx q[3];
rz(-1.937791) q[3];
sx q[3];
rz(1.7758241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0492964) q[2];
sx q[2];
rz(-0.80731097) q[2];
sx q[2];
rz(2.2558007) q[2];
rz(2.2057335) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7037999) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.6920775) q[0];
rz(1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(1.6922916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3291849) q[0];
sx q[0];
rz(-1.158876) q[0];
sx q[0];
rz(0.054188577) q[0];
x q[1];
rz(2.7141391) q[2];
sx q[2];
rz(-1.1910254) q[2];
sx q[2];
rz(0.7537656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76043452) q[1];
sx q[1];
rz(-1.0236003) q[1];
sx q[1];
rz(0.24914279) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0725771) q[3];
sx q[3];
rz(-0.3007362) q[3];
sx q[3];
rz(-1.6274483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8861822) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(-0.6748684) q[2];
rz(-2.1700962) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424727) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(-2.9091861) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(2.7276298) q[2];
sx q[2];
rz(-2.064075) q[2];
sx q[2];
rz(1.0309564) q[2];
rz(1.483199) q[3];
sx q[3];
rz(-2.4110553) q[3];
sx q[3];
rz(0.72881107) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
