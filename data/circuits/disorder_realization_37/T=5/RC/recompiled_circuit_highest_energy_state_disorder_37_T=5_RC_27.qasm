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
rz(2.2651894) q[0];
sx q[0];
rz(9.6776008) q[0];
rz(-1.9934935) q[1];
sx q[1];
rz(-0.91528457) q[1];
sx q[1];
rz(0.80438703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0641805) q[0];
sx q[0];
rz(-1.8837116) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43008606) q[2];
sx q[2];
rz(-1.7835435) q[2];
sx q[2];
rz(1.2141922) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8411257) q[1];
sx q[1];
rz(-1.2004832) q[1];
sx q[1];
rz(-1.200586) q[1];
rz(-pi) q[2];
rz(0.31320235) q[3];
sx q[3];
rz(-0.36010183) q[3];
sx q[3];
rz(1.859377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.09314166) q[2];
sx q[2];
rz(-2.1790049) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63808477) q[0];
sx q[0];
rz(-2.8758949) q[0];
sx q[0];
rz(-2.9124394) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1198657) q[0];
sx q[0];
rz(-1.3029126) q[0];
sx q[0];
rz(0.88587205) q[0];
rz(-pi) q[1];
rz(-0.83723703) q[2];
sx q[2];
rz(-1.9847437) q[2];
sx q[2];
rz(-2.4606103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9728127) q[1];
sx q[1];
rz(-0.75835627) q[1];
sx q[1];
rz(0.43557628) q[1];
x q[2];
rz(0.76477625) q[3];
sx q[3];
rz(-1.4524609) q[3];
sx q[3];
rz(1.8272994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68636346) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(2.6658106) q[2];
rz(2.0717715) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0723202) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(1.9092165) q[0];
rz(-0.45066372) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(0.88952363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43972507) q[0];
sx q[0];
rz(-0.98977463) q[0];
sx q[0];
rz(1.0807316) q[0];
rz(-2.1289775) q[2];
sx q[2];
rz(-2.0348437) q[2];
sx q[2];
rz(1.2166263) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3308308) q[1];
sx q[1];
rz(-2.2137801) q[1];
sx q[1];
rz(0.26442702) q[1];
rz(-pi) q[2];
rz(-2.7660288) q[3];
sx q[3];
rz(-1.2455539) q[3];
sx q[3];
rz(-0.90876814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-2.7165518) q[2];
sx q[2];
rz(1.8827776) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71037978) q[0];
sx q[0];
rz(-0.90573913) q[0];
sx q[0];
rz(1.6718965) q[0];
rz(-1.9944893) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(0.9009487) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6916136) q[0];
sx q[0];
rz(-1.2581345) q[0];
sx q[0];
rz(1.481196) q[0];
rz(0.6154279) q[2];
sx q[2];
rz(-2.0912898) q[2];
sx q[2];
rz(-2.8857638) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4576396) q[1];
sx q[1];
rz(-1.4767273) q[1];
sx q[1];
rz(-0.72503083) q[1];
x q[2];
rz(-2.9649686) q[3];
sx q[3];
rz(-2.400853) q[3];
sx q[3];
rz(-0.82308849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6441696) q[2];
sx q[2];
rz(-1.7449417) q[2];
sx q[2];
rz(-0.72745848) q[2];
rz(-0.71508956) q[3];
sx q[3];
rz(-0.46052027) q[3];
sx q[3];
rz(-2.6434744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(2.4053307) q[0];
rz(2.5212506) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(-1.5249407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10337457) q[0];
sx q[0];
rz(-0.056110121) q[0];
sx q[0];
rz(2.9655064) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9061766) q[2];
sx q[2];
rz(-1.7896626) q[2];
sx q[2];
rz(2.9218505) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1279631) q[1];
sx q[1];
rz(-1.1971601) q[1];
sx q[1];
rz(-2.29261) q[1];
rz(-pi) q[2];
rz(2.3518538) q[3];
sx q[3];
rz(-2.2135255) q[3];
sx q[3];
rz(1.5585043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6058495) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(-0.95174754) q[2];
rz(2.7239299) q[3];
sx q[3];
rz(-0.2499191) q[3];
sx q[3];
rz(1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34053892) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(0.75403768) q[0];
rz(0.89318371) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(2.975614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5170068) q[0];
sx q[0];
rz(-0.69147325) q[0];
sx q[0];
rz(-0.16384478) q[0];
rz(-pi) q[1];
rz(1.3909666) q[2];
sx q[2];
rz(-1.0713312) q[2];
sx q[2];
rz(1.7714372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80604078) q[1];
sx q[1];
rz(-1.1245831) q[1];
sx q[1];
rz(-0.70202804) q[1];
x q[2];
rz(-3.0021871) q[3];
sx q[3];
rz(-1.8434815) q[3];
sx q[3];
rz(0.26095873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(0.005391187) q[2];
rz(0.088717669) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(0.79571342) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8832815) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(2.2329277) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(3.0970792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712898) q[0];
sx q[0];
rz(-1.7006386) q[0];
sx q[0];
rz(-1.4128039) q[0];
rz(-0.037864121) q[2];
sx q[2];
rz(-1.711039) q[2];
sx q[2];
rz(-0.87272206) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14340885) q[1];
sx q[1];
rz(-2.1381731) q[1];
sx q[1];
rz(-0.25556775) q[1];
x q[2];
rz(-1.4511075) q[3];
sx q[3];
rz(-1.5482248) q[3];
sx q[3];
rz(-0.97609568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8160416) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(2.6204056) q[2];
rz(-0.19736396) q[3];
sx q[3];
rz(-1.7557996) q[3];
sx q[3];
rz(0.50281966) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40081438) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(-2.8863696) q[0];
rz(2.9995645) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(2.7878888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3981374) q[0];
sx q[0];
rz(-0.48851911) q[0];
sx q[0];
rz(2.3228881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92112215) q[2];
sx q[2];
rz(-0.73474681) q[2];
sx q[2];
rz(-0.91780797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28359828) q[1];
sx q[1];
rz(-1.9093175) q[1];
sx q[1];
rz(1.1095065) q[1];
x q[2];
rz(-0.74905101) q[3];
sx q[3];
rz(-2.0045677) q[3];
sx q[3];
rz(-1.6925136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8346617) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(1.5241415) q[2];
rz(-0.8664242) q[3];
sx q[3];
rz(-1.4640936) q[3];
sx q[3];
rz(2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50188142) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(-1.9068708) q[1];
sx q[1];
rz(-0.23209485) q[1];
sx q[1];
rz(0.51188767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8361266) q[0];
sx q[0];
rz(-2.9752467) q[0];
sx q[0];
rz(-1.2438828) q[0];
rz(-pi) q[1];
rz(-0.16361605) q[2];
sx q[2];
rz(-0.77218845) q[2];
sx q[2];
rz(2.8761187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3583957) q[1];
sx q[1];
rz(-1.8640717) q[1];
sx q[1];
rz(-0.44289987) q[1];
rz(2.7861106) q[3];
sx q[3];
rz(-1.0159229) q[3];
sx q[3];
rz(-1.9112089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(0.39719886) q[2];
rz(-1.0154137) q[3];
sx q[3];
rz(-1.0011324) q[3];
sx q[3];
rz(2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7427202) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(2.7040828) q[0];
rz(2.4028026) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(1.662426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6300059) q[0];
sx q[0];
rz(-1.6264722) q[0];
sx q[0];
rz(-0.48752174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.871256) q[2];
sx q[2];
rz(-2.766937) q[2];
sx q[2];
rz(-2.7596291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2228415) q[1];
sx q[1];
rz(-1.4645618) q[1];
sx q[1];
rz(1.5083471) q[1];
rz(2.4303834) q[3];
sx q[3];
rz(-2.2235302) q[3];
sx q[3];
rz(2.0026596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2207569) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-0.1753359) q[2];
rz(-1.0026503) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6709082) q[0];
sx q[0];
rz(-1.6841472) q[0];
sx q[0];
rz(2.0367784) q[0];
rz(2.9910174) q[1];
sx q[1];
rz(-2.2655948) q[1];
sx q[1];
rz(1.2949952) q[1];
rz(0.19828037) q[2];
sx q[2];
rz(-1.56326) q[2];
sx q[2];
rz(2.9247912) q[2];
rz(2.8908875) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
