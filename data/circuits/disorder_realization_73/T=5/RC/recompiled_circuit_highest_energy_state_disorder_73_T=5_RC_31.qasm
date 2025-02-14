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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(0.90176982) q[0];
rz(-1.8094485) q[1];
sx q[1];
rz(-0.44337115) q[1];
sx q[1];
rz(-2.3763357) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741824) q[0];
sx q[0];
rz(-1.2077792) q[0];
sx q[0];
rz(0.98591759) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95633219) q[2];
sx q[2];
rz(-2.7354089) q[2];
sx q[2];
rz(1.9910016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74435788) q[1];
sx q[1];
rz(-1.9089948) q[1];
sx q[1];
rz(1.6133135) q[1];
x q[2];
rz(1.2918858) q[3];
sx q[3];
rz(-1.4456914) q[3];
sx q[3];
rz(0.58663128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8454933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(-1.6695401) q[2];
rz(-1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(0.53372395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.958309) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(2.2928152) q[0];
rz(-3.1065885) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(-2.1880207) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6153107) q[0];
sx q[0];
rz(-2.1325975) q[0];
sx q[0];
rz(2.1192106) q[0];
rz(1.1814647) q[2];
sx q[2];
rz(-1.8539394) q[2];
sx q[2];
rz(-2.2623537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9417022) q[1];
sx q[1];
rz(-2.7459956) q[1];
sx q[1];
rz(2.0989749) q[1];
x q[2];
rz(-0.75462975) q[3];
sx q[3];
rz(-2.0622232) q[3];
sx q[3];
rz(1.400465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1635052) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(-0.84954849) q[2];
rz(-0.16544011) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(0.74976841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(-2.7396696) q[0];
rz(-0.15332128) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(1.1361928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83602609) q[0];
sx q[0];
rz(-1.5587806) q[0];
sx q[0];
rz(1.9703321) q[0];
x q[1];
rz(-0.48353124) q[2];
sx q[2];
rz(-1.8147959) q[2];
sx q[2];
rz(2.2738843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2991645) q[1];
sx q[1];
rz(-1.8106726) q[1];
sx q[1];
rz(-0.23934083) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4506684) q[3];
sx q[3];
rz(-2.5657585) q[3];
sx q[3];
rz(-1.3802647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0574657) q[2];
sx q[2];
rz(-2.2630313) q[2];
sx q[2];
rz(0.055056661) q[2];
rz(1.7940686) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-2.1373035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771626) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5675911) q[0];
rz(-2.4500627) q[1];
sx q[1];
rz(-0.66929308) q[1];
sx q[1];
rz(0.18542586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90339336) q[0];
sx q[0];
rz(-1.967634) q[0];
sx q[0];
rz(-1.7522041) q[0];
x q[1];
rz(-1.144997) q[2];
sx q[2];
rz(-1.123872) q[2];
sx q[2];
rz(1.9831234) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37764964) q[1];
sx q[1];
rz(-0.40622207) q[1];
sx q[1];
rz(2.2561314) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3641649) q[3];
sx q[3];
rz(-2.1071599) q[3];
sx q[3];
rz(0.11851507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4039679) q[2];
sx q[2];
rz(-1.7125968) q[2];
sx q[2];
rz(0.0052304012) q[2];
rz(-0.38749203) q[3];
sx q[3];
rz(-2.6369075) q[3];
sx q[3];
rz(-1.3030049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72652793) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(2.5256185) q[0];
rz(-1.9527324) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(0.51330769) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270638) q[0];
sx q[0];
rz(-1.5634057) q[0];
sx q[0];
rz(1.8523995) q[0];
x q[1];
rz(-0.31536021) q[2];
sx q[2];
rz(-0.7509481) q[2];
sx q[2];
rz(-0.89799228) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23780248) q[1];
sx q[1];
rz(-0.93751838) q[1];
sx q[1];
rz(1.4813745) q[1];
rz(-pi) q[2];
rz(-1.1034455) q[3];
sx q[3];
rz(-1.4650405) q[3];
sx q[3];
rz(2.8018564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1335699) q[2];
sx q[2];
rz(-1.8654537) q[2];
sx q[2];
rz(-0.7424773) q[2];
rz(-1.4437458) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(2.0402563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3105069) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(0.35980862) q[0];
rz(0.85044914) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(0.46583072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208) q[0];
sx q[0];
rz(-0.75220901) q[0];
sx q[0];
rz(-1.2515995) q[0];
rz(1.7066782) q[2];
sx q[2];
rz(-2.3294221) q[2];
sx q[2];
rz(3.0480609) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40738525) q[1];
sx q[1];
rz(-1.3300597) q[1];
sx q[1];
rz(-1.5182785) q[1];
x q[2];
rz(0.96420607) q[3];
sx q[3];
rz(-1.0830942) q[3];
sx q[3];
rz(1.4473947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.74238527) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(-2.7772389) q[2];
rz(-0.24389167) q[3];
sx q[3];
rz(-3.0920005) q[3];
sx q[3];
rz(-1.9243139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(1.97557) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(-2.7551415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7313663) q[0];
sx q[0];
rz(-1.7803734) q[0];
sx q[0];
rz(-3.1102577) q[0];
rz(-pi) q[1];
rz(-2.4555444) q[2];
sx q[2];
rz(-0.39271388) q[2];
sx q[2];
rz(-0.60143747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7658479) q[1];
sx q[1];
rz(-1.6465636) q[1];
sx q[1];
rz(0.025546207) q[1];
rz(1.7771882) q[3];
sx q[3];
rz(-0.55183059) q[3];
sx q[3];
rz(-1.3760374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7361136) q[2];
sx q[2];
rz(-1.6454641) q[2];
sx q[2];
rz(-0.22769895) q[2];
rz(-2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(-2.8480215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-2.7135799) q[0];
rz(2.1233066) q[1];
sx q[1];
rz(-0.4363474) q[1];
sx q[1];
rz(-2.2705618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7334325) q[0];
sx q[0];
rz(-2.1299612) q[0];
sx q[0];
rz(-0.59691043) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0592878) q[2];
sx q[2];
rz(-2.7370484) q[2];
sx q[2];
rz(1.900857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5586413) q[1];
sx q[1];
rz(-1.691733) q[1];
sx q[1];
rz(-1.7128829) q[1];
rz(-pi) q[2];
rz(1.0139129) q[3];
sx q[3];
rz(-0.96864163) q[3];
sx q[3];
rz(-2.7891261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13228664) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(1.8629249) q[2];
rz(-0.57847413) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(1.1540029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4671675) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(-0.59980741) q[0];
rz(1.2990052) q[1];
sx q[1];
rz(-1.6103585) q[1];
sx q[1];
rz(-0.64186796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917832) q[0];
sx q[0];
rz(-1.9115752) q[0];
sx q[0];
rz(-2.6343915) q[0];
rz(-pi) q[1];
rz(-2.3360905) q[2];
sx q[2];
rz(-1.6751211) q[2];
sx q[2];
rz(0.25539216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.95002) q[1];
sx q[1];
rz(-2.2784462) q[1];
sx q[1];
rz(-1.5465082) q[1];
rz(0.55223744) q[3];
sx q[3];
rz(-0.83697666) q[3];
sx q[3];
rz(-0.53572922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.32144) q[2];
sx q[2];
rz(-1.7165311) q[2];
sx q[2];
rz(0.610262) q[2];
rz(-0.52014703) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(2.648073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34058061) q[0];
sx q[0];
rz(-1.8980674) q[0];
sx q[0];
rz(-2.2035759) q[0];
rz(0.34171379) q[1];
sx q[1];
rz(-1.382117) q[1];
sx q[1];
rz(2.9843073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8167687) q[0];
sx q[0];
rz(-0.95476645) q[0];
sx q[0];
rz(1.8094713) q[0];
rz(0.089265169) q[2];
sx q[2];
rz(-2.3815739) q[2];
sx q[2];
rz(-1.6940728) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57676901) q[1];
sx q[1];
rz(-1.8722459) q[1];
sx q[1];
rz(-0.5524854) q[1];
rz(0.25453849) q[3];
sx q[3];
rz(-1.1738281) q[3];
sx q[3];
rz(-2.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9226795) q[2];
sx q[2];
rz(-0.34505406) q[2];
sx q[2];
rz(1.2738796) q[2];
rz(-3.1359361) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(-0.98102942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2904749) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(-0.46650096) q[1];
sx q[1];
rz(-1.7955045) q[1];
sx q[1];
rz(-2.4972965) q[1];
rz(-0.90032719) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(0.43177035) q[3];
sx q[3];
rz(-1.645586) q[3];
sx q[3];
rz(0.077787568) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
