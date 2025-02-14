OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(-2.9105817) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92316931) q[0];
sx q[0];
rz(-1.8418663) q[0];
sx q[0];
rz(2.0725033) q[0];
rz(-1.3363935) q[2];
sx q[2];
rz(-1.3318828) q[2];
sx q[2];
rz(-0.51942458) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8637064) q[1];
sx q[1];
rz(-1.7369629) q[1];
sx q[1];
rz(-1.4708797) q[1];
x q[2];
rz(-1.1679959) q[3];
sx q[3];
rz(-0.19013817) q[3];
sx q[3];
rz(-2.8791752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7370558) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-0.36049584) q[2];
rz(-1.1245842) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8812113) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(2.4847109) q[0];
rz(0.33292133) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(0.93516707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012197709) q[0];
sx q[0];
rz(-3.0327065) q[0];
sx q[0];
rz(2.8588141) q[0];
rz(-pi) q[1];
rz(2.8843811) q[2];
sx q[2];
rz(-1.6597444) q[2];
sx q[2];
rz(2.6572029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6015687) q[1];
sx q[1];
rz(-1.4238796) q[1];
sx q[1];
rz(-1.6107035) q[1];
rz(-pi) q[2];
rz(2.786119) q[3];
sx q[3];
rz(-1.7661816) q[3];
sx q[3];
rz(0.030578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3071345) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(-1.9251941) q[2];
rz(-2.6136716) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6127748) q[0];
sx q[0];
rz(-0.82614326) q[0];
sx q[0];
rz(-2.5566027) q[0];
rz(-3.0168369) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(-1.4488719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.751469) q[0];
sx q[0];
rz(-1.5737783) q[0];
sx q[0];
rz(-3.1389159) q[0];
rz(-pi) q[1];
rz(3.0694088) q[2];
sx q[2];
rz(-1.584314) q[2];
sx q[2];
rz(1.1987606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6022012) q[1];
sx q[1];
rz(-0.81967205) q[1];
sx q[1];
rz(-0.96266268) q[1];
rz(2.8040941) q[3];
sx q[3];
rz(-1.322116) q[3];
sx q[3];
rz(1.6832222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8831732) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(1.4460571) q[2];
rz(0.7575194) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739173) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(2.0539334) q[0];
rz(2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-2.1813724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86726928) q[0];
sx q[0];
rz(-1.3057297) q[0];
sx q[0];
rz(0.041170711) q[0];
rz(-1.6702508) q[2];
sx q[2];
rz(-1.5410454) q[2];
sx q[2];
rz(-3.0029675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31611495) q[1];
sx q[1];
rz(-1.1255472) q[1];
sx q[1];
rz(-0.47834217) q[1];
rz(0.26953617) q[3];
sx q[3];
rz(-1.4785267) q[3];
sx q[3];
rz(0.59449457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20597657) q[2];
sx q[2];
rz(-1.874186) q[2];
sx q[2];
rz(-1.6440294) q[2];
rz(-2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(-2.5373051) q[0];
rz(-1.1445649) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(-0.42246517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9086994) q[0];
sx q[0];
rz(-0.18747231) q[0];
sx q[0];
rz(2.0220246) q[0];
x q[1];
rz(2.3461012) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(2.7249641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8026655) q[1];
sx q[1];
rz(-0.69310729) q[1];
sx q[1];
rz(-3.0285804) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72679584) q[3];
sx q[3];
rz(-0.49203792) q[3];
sx q[3];
rz(-0.43471042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8233238) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(-2.4853415) q[2];
rz(-1.6107791) q[3];
sx q[3];
rz(-1.8717513) q[3];
sx q[3];
rz(-0.58073616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-2.1298213) q[0];
sx q[0];
rz(-2.5166125) q[0];
rz(2.3896353) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(1.6065074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84941712) q[0];
sx q[0];
rz(-1.2818953) q[0];
sx q[0];
rz(1.6221415) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1691323) q[2];
sx q[2];
rz(-2.1876799) q[2];
sx q[2];
rz(-1.0765558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5175486) q[1];
sx q[1];
rz(-2.4772948) q[1];
sx q[1];
rz(2.7655927) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1702437) q[3];
sx q[3];
rz(-1.1068212) q[3];
sx q[3];
rz(-2.6148207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9525166) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(-0.9291741) q[2];
rz(-2.3164228) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(-1.1024124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(-1.7671385) q[0];
rz(-0.51042405) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-2.5183831) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7316298) q[0];
sx q[0];
rz(-1.4193168) q[0];
sx q[0];
rz(1.9069457) q[0];
x q[1];
rz(1.7665777) q[2];
sx q[2];
rz(-2.0718241) q[2];
sx q[2];
rz(-2.632338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53155073) q[1];
sx q[1];
rz(-1.4248214) q[1];
sx q[1];
rz(3.0290305) q[1];
rz(-pi) q[2];
rz(1.3680063) q[3];
sx q[3];
rz(-0.8335127) q[3];
sx q[3];
rz(1.1356419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0509384) q[2];
sx q[2];
rz(-2.3263558) q[2];
sx q[2];
rz(1.1673002) q[2];
rz(0.022631571) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24564329) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(-2.1242712) q[0];
rz(-1.7474878) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(-2.5140433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.548259) q[0];
sx q[0];
rz(-0.98932668) q[0];
sx q[0];
rz(2.6618746) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6553788) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(1.6135474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4269451) q[1];
sx q[1];
rz(-1.3332187) q[1];
sx q[1];
rz(1.7725043) q[1];
x q[2];
rz(-1.2305607) q[3];
sx q[3];
rz(-2.0215394) q[3];
sx q[3];
rz(-2.4204202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5158186) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-2.7308357) q[2];
rz(-3.0913894) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5589767) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(-0.26891747) q[0];
rz(-0.31002632) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(-0.1300098) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63142473) q[0];
sx q[0];
rz(-1.893259) q[0];
sx q[0];
rz(2.9084951) q[0];
rz(-pi) q[1];
rz(0.35138826) q[2];
sx q[2];
rz(-1.2341502) q[2];
sx q[2];
rz(1.2996246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42327729) q[1];
sx q[1];
rz(-2.1287103) q[1];
sx q[1];
rz(1.3807135) q[1];
x q[2];
rz(2.1653529) q[3];
sx q[3];
rz(-2.3701982) q[3];
sx q[3];
rz(0.4679799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2076063) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(0.096573528) q[2];
rz(-1.5038331) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132975) q[0];
sx q[0];
rz(-2.9533563) q[0];
sx q[0];
rz(2.6085594) q[0];
rz(-3.10532) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(1.4520377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55646) q[0];
sx q[0];
rz(-2.1341354) q[0];
sx q[0];
rz(-0.4507395) q[0];
rz(-pi) q[1];
rz(-1.6797941) q[2];
sx q[2];
rz(-1.8882635) q[2];
sx q[2];
rz(0.017680971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.888537) q[1];
sx q[1];
rz(-1.5408684) q[1];
sx q[1];
rz(2.7821343) q[1];
rz(-pi) q[2];
rz(1.6598654) q[3];
sx q[3];
rz(-2.0348106) q[3];
sx q[3];
rz(-2.9651412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4286246) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(-0.0326322) q[2];
rz(-0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(-2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(-2.3445917) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-0.47427872) q[2];
sx q[2];
rz(-2.5534292) q[2];
sx q[2];
rz(-3.0234887) q[2];
rz(2.3305994) q[3];
sx q[3];
rz(-1.4832433) q[3];
sx q[3];
rz(1.087838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
