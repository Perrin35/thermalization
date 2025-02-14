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
rz(2.8887698) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(-2.2263081) q[1];
sx q[1];
rz(-0.80438703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0067449) q[0];
sx q[0];
rz(-2.6383556) q[0];
sx q[0];
rz(-2.2630967) q[0];
rz(0.43008606) q[2];
sx q[2];
rz(-1.3580492) q[2];
sx q[2];
rz(1.2141922) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6614395) q[1];
sx q[1];
rz(-2.6240908) q[1];
sx q[1];
rz(0.75019849) q[1];
x q[2];
rz(-0.34396307) q[3];
sx q[3];
rz(-1.6795782) q[3];
sx q[3];
rz(-0.0057084486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.09314166) q[2];
sx q[2];
rz(-0.96258771) q[2];
sx q[2];
rz(-2.2691881) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.395697) q[3];
sx q[3];
rz(-0.083757639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5035079) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(-0.22915325) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.3399905) q[1];
sx q[1];
rz(2.4280039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23799831) q[0];
sx q[0];
rz(-0.9146713) q[0];
sx q[0];
rz(0.34060069) q[0];
x q[1];
rz(-0.83723703) q[2];
sx q[2];
rz(-1.9847437) q[2];
sx q[2];
rz(0.6809823) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9728127) q[1];
sx q[1];
rz(-2.3832364) q[1];
sx q[1];
rz(-0.43557628) q[1];
rz(-pi) q[2];
rz(1.7341033) q[3];
sx q[3];
rz(-2.3288832) q[3];
sx q[3];
rz(0.14369609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-0.31193048) q[2];
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
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0723202) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(-0.45066372) q[1];
sx q[1];
rz(-1.9602937) q[1];
sx q[1];
rz(0.88952363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33143932) q[0];
sx q[0];
rz(-2.4001849) q[0];
sx q[0];
rz(-2.5196694) q[0];
rz(-pi) q[1];
rz(0.53307791) q[2];
sx q[2];
rz(-2.0642082) q[2];
sx q[2];
rz(-0.62668884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8107618) q[1];
sx q[1];
rz(-0.92781258) q[1];
sx q[1];
rz(2.8771656) q[1];
x q[2];
rz(2.7660288) q[3];
sx q[3];
rz(-1.8960388) q[3];
sx q[3];
rz(2.2328245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.687261) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(1.8827776) q[2];
rz(-1.9076294) q[3];
sx q[3];
rz(-2.0089269) q[3];
sx q[3];
rz(-3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312129) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(-1.4696962) q[0];
rz(-1.9944893) q[1];
sx q[1];
rz(-2.1397782) q[1];
sx q[1];
rz(0.9009487) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44997901) q[0];
sx q[0];
rz(-1.8834582) q[0];
sx q[0];
rz(1.6603966) q[0];
x q[1];
rz(0.9587165) q[2];
sx q[2];
rz(-1.0462073) q[2];
sx q[2];
rz(-2.1647349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.03012499) q[1];
sx q[1];
rz(-2.2919167) q[1];
sx q[1];
rz(-1.6961873) q[1];
x q[2];
rz(2.4086558) q[3];
sx q[3];
rz(-1.4519435) q[3];
sx q[3];
rz(0.87866966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(2.4141342) q[2];
rz(2.4265031) q[3];
sx q[3];
rz(-2.6810724) q[3];
sx q[3];
rz(-0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6775976) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(-0.736262) q[0];
rz(-2.5212506) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(-1.6166519) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10337457) q[0];
sx q[0];
rz(-0.056110121) q[0];
sx q[0];
rz(-2.9655064) q[0];
x q[1];
rz(-1.235416) q[2];
sx q[2];
rz(-1.35193) q[2];
sx q[2];
rz(-0.21974213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.013629524) q[1];
sx q[1];
rz(-1.1971601) q[1];
sx q[1];
rz(2.29261) q[1];
x q[2];
rz(2.3297263) q[3];
sx q[3];
rz(-2.1694739) q[3];
sx q[3];
rz(-2.6177399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6058495) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(0.95174754) q[2];
rz(2.7239299) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(-1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34053892) q[0];
sx q[0];
rz(-1.2507573) q[0];
sx q[0];
rz(0.75403768) q[0];
rz(-2.2484089) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(-0.16597861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9611325) q[0];
sx q[0];
rz(-1.6749973) q[0];
sx q[0];
rz(-2.4567338) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31678172) q[2];
sx q[2];
rz(-0.52826476) q[2];
sx q[2];
rz(2.1342017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7268035) q[1];
sx q[1];
rz(-0.94894743) q[1];
sx q[1];
rz(1.0110823) q[1];
rz(1.1096438) q[3];
sx q[3];
rz(-0.30546092) q[3];
sx q[3];
rz(2.4003254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7665427) q[2];
sx q[2];
rz(-1.1376209) q[2];
sx q[2];
rz(-0.005391187) q[2];
rz(3.052875) q[3];
sx q[3];
rz(-1.5327449) q[3];
sx q[3];
rz(-0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2583112) q[0];
sx q[0];
rz(-1.9323876) q[0];
sx q[0];
rz(-0.2051556) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7617247) q[0];
sx q[0];
rz(-1.7274478) q[0];
sx q[0];
rz(3.0101315) q[0];
rz(1.3088063) q[2];
sx q[2];
rz(-0.14523187) q[2];
sx q[2];
rz(2.0042208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5959357) q[1];
sx q[1];
rz(-0.6164729) q[1];
sx q[1];
rz(-1.9484387) q[1];
x q[2];
rz(0.022734108) q[3];
sx q[3];
rz(-1.6904545) q[3];
sx q[3];
rz(-2.5496063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32555106) q[2];
sx q[2];
rz(-2.6191235) q[2];
sx q[2];
rz(2.6204056) q[2];
rz(-0.19736396) q[3];
sx q[3];
rz(-1.7557996) q[3];
sx q[3];
rz(-2.638773) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40081438) q[0];
sx q[0];
rz(-2.9476808) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(-2.9995645) q[1];
sx q[1];
rz(-1.7250926) q[1];
sx q[1];
rz(-2.7878888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557705) q[0];
sx q[0];
rz(-1.9206128) q[0];
sx q[0];
rz(-2.7932998) q[0];
rz(2.2204705) q[2];
sx q[2];
rz(-2.4068458) q[2];
sx q[2];
rz(0.91780797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28359828) q[1];
sx q[1];
rz(-1.9093175) q[1];
sx q[1];
rz(-1.1095065) q[1];
rz(-pi) q[2];
rz(-2.5442637) q[3];
sx q[3];
rz(-0.84395614) q[3];
sx q[3];
rz(-2.8387808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3069309) q[2];
sx q[2];
rz(-0.31535172) q[2];
sx q[2];
rz(-1.5241415) q[2];
rz(2.2751685) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(0.58969897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50188142) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(-1.7891275) q[0];
rz(-1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-0.51188767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8361266) q[0];
sx q[0];
rz(-0.16634596) q[0];
sx q[0];
rz(1.2438828) q[0];
rz(1.7281248) q[2];
sx q[2];
rz(-2.3300588) q[2];
sx q[2];
rz(0.038977101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3583957) q[1];
sx q[1];
rz(-1.8640717) q[1];
sx q[1];
rz(2.6986928) q[1];
x q[2];
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
rz(-1.1341683) q[2];
sx q[2];
rz(-1.8210501) q[2];
sx q[2];
rz(-2.7443938) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
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
rz(-0.80016017) q[1];
sx q[1];
rz(1.4791666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0453607) q[0];
sx q[0];
rz(-0.49043617) q[0];
sx q[0];
rz(3.0231721) q[0];
rz(1.9301038) q[2];
sx q[2];
rz(-1.4622765) q[2];
sx q[2];
rz(1.6720275) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4829172) q[1];
sx q[1];
rz(-1.5086996) q[1];
sx q[1];
rz(-3.0351522) q[1];
rz(-pi) q[2];
rz(-0.71120925) q[3];
sx q[3];
rz(-0.91806245) q[3];
sx q[3];
rz(1.1389331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9208357) q[2];
sx q[2];
rz(-1.3584542) q[2];
sx q[2];
rz(-2.9662568) q[2];
rz(-2.1389424) q[3];
sx q[3];
rz(-0.23313871) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709082) q[0];
sx q[0];
rz(-1.4574454) q[0];
sx q[0];
rz(-1.1048143) q[0];
rz(0.15057527) q[1];
sx q[1];
rz(-0.87599788) q[1];
sx q[1];
rz(-1.8465975) q[1];
rz(2.9433123) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
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
