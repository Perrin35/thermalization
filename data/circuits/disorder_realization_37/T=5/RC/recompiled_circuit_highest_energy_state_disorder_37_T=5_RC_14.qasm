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
rz(-0.77684075) q[0];
sx q[0];
rz(-0.87640327) q[0];
sx q[0];
rz(-0.25282282) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(4.0568772) q[1];
sx q[1];
rz(8.6203909) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3768965) q[0];
sx q[0];
rz(-1.1904556) q[0];
sx q[0];
rz(2.8036731) q[0];
x q[1];
rz(-1.3374653) q[2];
sx q[2];
rz(-1.1510282) q[2];
sx q[2];
rz(0.45316089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.731718) q[1];
sx q[1];
rz(-1.2267641) q[1];
sx q[1];
rz(2.7469977) q[1];
x q[2];
rz(-1.4553045) q[3];
sx q[3];
rz(-1.9126429) q[3];
sx q[3];
rz(-1.6153743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.048451) q[2];
sx q[2];
rz(-2.1790049) q[2];
sx q[2];
rz(-2.2691881) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63808477) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(-2.9124394) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9035943) q[0];
sx q[0];
rz(-2.2269214) q[0];
sx q[0];
rz(2.800992) q[0];
rz(-pi) q[1];
rz(-2.607478) q[2];
sx q[2];
rz(-2.2306109) q[2];
sx q[2];
rz(0.5420064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9728127) q[1];
sx q[1];
rz(-0.75835627) q[1];
sx q[1];
rz(-0.43557628) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3768164) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-0.31193048) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(-1.0698211) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(-2.3750335) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0692724) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(-0.88952363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7018676) q[0];
sx q[0];
rz(-0.98977463) q[0];
sx q[0];
rz(-2.0608611) q[0];
x q[1];
rz(-2.6085147) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(-0.62668884) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3865859) q[1];
sx q[1];
rz(-0.68803794) q[1];
sx q[1];
rz(-1.2351456) q[1];
rz(-pi) q[2];
rz(-2.398185) q[3];
sx q[3];
rz(-2.6498859) q[3];
sx q[3];
rz(-1.3431026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4543317) q[2];
sx q[2];
rz(-2.7165518) q[2];
sx q[2];
rz(1.2588151) q[2];
rz(-1.2339633) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.71037978) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(-1.6718965) q[0];
rz(1.1471033) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(0.9009487) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44997901) q[0];
sx q[0];
rz(-1.2581345) q[0];
sx q[0];
rz(-1.481196) q[0];
rz(-pi) q[1];
rz(0.9587165) q[2];
sx q[2];
rz(-2.0953853) q[2];
sx q[2];
rz(2.1647349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1114677) q[1];
sx q[1];
rz(-2.2919167) q[1];
sx q[1];
rz(-1.6961873) q[1];
rz(0.73293682) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(2.262923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6441696) q[2];
sx q[2];
rz(-1.7449417) q[2];
sx q[2];
rz(-0.72745848) q[2];
rz(-2.4265031) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46399507) q[0];
sx q[0];
rz(-1.7514739) q[0];
sx q[0];
rz(-0.736262) q[0];
rz(-0.62034208) q[1];
sx q[1];
rz(-2.5963929) q[1];
sx q[1];
rz(-1.6166519) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0382181) q[0];
sx q[0];
rz(-0.056110121) q[0];
sx q[0];
rz(-2.9655064) q[0];
rz(2.1651046) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(1.9081685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1641008) q[1];
sx q[1];
rz(-0.79704801) q[1];
sx q[1];
rz(-2.1063095) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81186632) q[3];
sx q[3];
rz(-0.97211876) q[3];
sx q[3];
rz(2.6177399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5357431) q[2];
sx q[2];
rz(-0.77631408) q[2];
sx q[2];
rz(0.95174754) q[2];
rz(-2.7239299) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(-1.9865659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053892) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(2.387555) q[0];
rz(-2.2484089) q[1];
sx q[1];
rz(-2.1008396) q[1];
sx q[1];
rz(0.16597861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5170068) q[0];
sx q[0];
rz(-0.69147325) q[0];
sx q[0];
rz(-0.16384478) q[0];
x q[1];
rz(-1.3909666) q[2];
sx q[2];
rz(-2.0702614) q[2];
sx q[2];
rz(1.7714372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80604078) q[1];
sx q[1];
rz(-1.1245831) q[1];
sx q[1];
rz(2.4395646) q[1];
rz(-pi) q[2];
rz(1.8460205) q[3];
sx q[3];
rz(-1.7050192) q[3];
sx q[3];
rz(1.272066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7665427) q[2];
sx q[2];
rz(-1.1376209) q[2];
sx q[2];
rz(3.1362015) q[2];
rz(-0.088717669) q[3];
sx q[3];
rz(-1.6088477) q[3];
sx q[3];
rz(0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8832815) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(2.9364371) q[0];
rz(-2.2329277) q[1];
sx q[1];
rz(-0.87037194) q[1];
sx q[1];
rz(3.0970792) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7617247) q[0];
sx q[0];
rz(-1.4141448) q[0];
sx q[0];
rz(-3.0101315) q[0];
rz(-1.3088063) q[2];
sx q[2];
rz(-2.9963608) q[2];
sx q[2];
rz(-1.1373718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8537112) q[1];
sx q[1];
rz(-1.3559623) q[1];
sx q[1];
rz(0.98835215) q[1];
rz(1.4511075) q[3];
sx q[3];
rz(-1.5482248) q[3];
sx q[3];
rz(-2.165497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8160416) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(-0.52118707) q[2];
rz(0.19736396) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-2.638773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.40081438) q[0];
sx q[0];
rz(-2.9476808) q[0];
sx q[0];
rz(-2.8863696) q[0];
rz(-2.9995645) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(-2.7878888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784894) q[0];
sx q[0];
rz(-1.8971839) q[0];
sx q[0];
rz(1.9410237) q[0];
x q[1];
rz(2.1944517) q[2];
sx q[2];
rz(-1.9883687) q[2];
sx q[2];
rz(3.0018978) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6981909) q[1];
sx q[1];
rz(-2.576722) q[1];
sx q[1];
rz(-2.2400212) q[1];
x q[2];
rz(2.1347624) q[3];
sx q[3];
rz(-0.90475291) q[3];
sx q[3];
rz(-2.6473235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3069309) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(-1.6174512) q[2];
rz(0.8664242) q[3];
sx q[3];
rz(-1.4640936) q[3];
sx q[3];
rz(-2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6397112) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(1.3524652) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-2.629705) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0842132) q[0];
sx q[0];
rz(-1.6239927) q[0];
sx q[0];
rz(-1.7284815) q[0];
x q[1];
rz(2.9779766) q[2];
sx q[2];
rz(-2.3694042) q[2];
sx q[2];
rz(0.26547394) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7929188) q[1];
sx q[1];
rz(-1.9935384) q[1];
sx q[1];
rz(-1.2482367) q[1];
rz(2.0824349) q[3];
sx q[3];
rz(-0.64877227) q[3];
sx q[3];
rz(1.8442475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1341683) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(0.39719886) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(0.5526244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7427202) q[0];
sx q[0];
rz(-1.2074559) q[0];
sx q[0];
rz(2.7040828) q[0];
rz(-0.73879009) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(1.662426) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2302983) q[0];
sx q[0];
rz(-1.084096) q[0];
sx q[0];
rz(-1.6337956) q[0];
rz(-pi) q[1];
rz(-3.0257341) q[2];
sx q[2];
rz(-1.9278952) q[2];
sx q[2];
rz(0.060572123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2228415) q[1];
sx q[1];
rz(-1.4645618) q[1];
sx q[1];
rz(-1.5083471) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4303834) q[3];
sx q[3];
rz(-2.2235302) q[3];
sx q[3];
rz(1.1389331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2207569) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-0.1753359) q[2];
rz(1.0026503) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.6709082) q[0];
sx q[0];
rz(-1.6841472) q[0];
sx q[0];
rz(2.0367784) q[0];
rz(-0.15057527) q[1];
sx q[1];
rz(-2.2655948) q[1];
sx q[1];
rz(1.2949952) q[1];
rz(3.1033517) q[2];
sx q[2];
rz(-2.943171) q[2];
sx q[2];
rz(1.3914862) q[2];
rz(-2.3948432) q[3];
sx q[3];
rz(-0.33573739) q[3];
sx q[3];
rz(2.4949069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
