OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(-2.28595) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(0.56245437) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320975) q[0];
sx q[0];
rz(-1.9350855) q[0];
sx q[0];
rz(3.0845736) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0777656) q[2];
sx q[2];
rz(-1.7071144) q[2];
sx q[2];
rz(2.9796114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89002883) q[1];
sx q[1];
rz(-1.9891796) q[1];
sx q[1];
rz(-2.0189925) q[1];
rz(-2.7482618) q[3];
sx q[3];
rz(-1.1926023) q[3];
sx q[3];
rz(-2.8761169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(2.3423024) q[2];
rz(2.6702787) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(2.2057064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039463194) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.348173) q[0];
rz(1.3735636) q[1];
sx q[1];
rz(-1.9919688) q[1];
sx q[1];
rz(1.1522393) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2179759) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(-0.38461916) q[0];
rz(2.466408) q[2];
sx q[2];
rz(-1.0567046) q[2];
sx q[2];
rz(-0.49066431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5785529) q[1];
sx q[1];
rz(-2.4336947) q[1];
sx q[1];
rz(-2.8824174) q[1];
rz(-pi) q[2];
rz(-2.3186734) q[3];
sx q[3];
rz(-2.4429818) q[3];
sx q[3];
rz(2.5320092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(-3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057864144) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-2.6237543) q[0];
rz(-0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(-2.6944366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0797821) q[0];
sx q[0];
rz(-1.634565) q[0];
sx q[0];
rz(-0.26173862) q[0];
rz(-pi) q[1];
rz(-2.4141623) q[2];
sx q[2];
rz(-2.0038249) q[2];
sx q[2];
rz(-2.7114781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4616926) q[1];
sx q[1];
rz(-2.5352806) q[1];
sx q[1];
rz(1.8099422) q[1];
rz(1.2918043) q[3];
sx q[3];
rz(-2.0800262) q[3];
sx q[3];
rz(-1.4470456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7670224) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(0.53263295) q[2];
rz(-2.8940708) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240876) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(-0.10661539) q[0];
rz(-2.7952349) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(1.9812298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713523) q[0];
sx q[0];
rz(-1.3246857) q[0];
sx q[0];
rz(-1.3918124) q[0];
rz(-0.71765064) q[2];
sx q[2];
rz(-2.1776878) q[2];
sx q[2];
rz(0.90536149) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6939099) q[1];
sx q[1];
rz(-2.6676148) q[1];
sx q[1];
rz(-0.86556566) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87507803) q[3];
sx q[3];
rz(-0.64662537) q[3];
sx q[3];
rz(0.88133206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0059263) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(2.5557319) q[3];
sx q[3];
rz(-1.9867089) q[3];
sx q[3];
rz(1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0214486) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(2.7914877) q[0];
rz(-2.1856951) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(-1.2299445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1317539) q[0];
sx q[0];
rz(-1.810605) q[0];
sx q[0];
rz(1.9480223) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7361264) q[2];
sx q[2];
rz(-1.606719) q[2];
sx q[2];
rz(-1.3816116) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38795162) q[1];
sx q[1];
rz(-1.759583) q[1];
sx q[1];
rz(1.2354047) q[1];
rz(-pi) q[2];
rz(2.941972) q[3];
sx q[3];
rz(-0.99310447) q[3];
sx q[3];
rz(-2.9132089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2436287) q[2];
sx q[2];
rz(-0.25883365) q[2];
sx q[2];
rz(-1.6492856) q[2];
rz(0.40488511) q[3];
sx q[3];
rz(-0.76571524) q[3];
sx q[3];
rz(0.83612061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355857) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(2.5591922) q[0];
rz(0.68663418) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.2581717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4674038) q[0];
sx q[0];
rz(-1.6961401) q[0];
sx q[0];
rz(2.0095429) q[0];
rz(-1.0354543) q[2];
sx q[2];
rz(-1.8796433) q[2];
sx q[2];
rz(0.18986407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2542779) q[1];
sx q[1];
rz(-1.6266903) q[1];
sx q[1];
rz(-0.069737597) q[1];
rz(-2.9223676) q[3];
sx q[3];
rz(-0.98837438) q[3];
sx q[3];
rz(-1.5150013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3758292) q[2];
sx q[2];
rz(-1.148369) q[2];
sx q[2];
rz(2.1235535) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.3956416) q[3];
sx q[3];
rz(1.6233981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(0.14353453) q[1];
sx q[1];
rz(-0.48964557) q[1];
sx q[1];
rz(0.20763436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26941368) q[0];
sx q[0];
rz(-2.2996443) q[0];
sx q[0];
rz(-0.4704041) q[0];
rz(-pi) q[1];
rz(2.1499499) q[2];
sx q[2];
rz(-1.8657547) q[2];
sx q[2];
rz(2.7911012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5112111) q[1];
sx q[1];
rz(-2.2464754) q[1];
sx q[1];
rz(2.8065794) q[1];
x q[2];
rz(-2.7782441) q[3];
sx q[3];
rz(-1.4814113) q[3];
sx q[3];
rz(-0.18401981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(2.5353954) q[2];
rz(2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(-0.70639759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48924482) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(1.0580753) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(-1.441997) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88272754) q[0];
sx q[0];
rz(-1.2611212) q[0];
sx q[0];
rz(-0.41559269) q[0];
rz(1.2953561) q[2];
sx q[2];
rz(-2.0552962) q[2];
sx q[2];
rz(-0.4415919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7322526) q[1];
sx q[1];
rz(-2.2870018) q[1];
sx q[1];
rz(0.90998896) q[1];
x q[2];
rz(-2.2820246) q[3];
sx q[3];
rz(-2.7894434) q[3];
sx q[3];
rz(-1.5052049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66120061) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(-1.2083758) q[2];
rz(0.66323534) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(1.2164345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709568) q[0];
sx q[0];
rz(-2.1747776) q[0];
sx q[0];
rz(3.048625) q[0];
rz(-1.2804821) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(-3.0063937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7057892) q[0];
sx q[0];
rz(-0.069878526) q[0];
sx q[0];
rz(2.4326434) q[0];
rz(-pi) q[1];
rz(-2.8835758) q[2];
sx q[2];
rz(-1.4224367) q[2];
sx q[2];
rz(0.49000636) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7635423) q[1];
sx q[1];
rz(-1.5494746) q[1];
sx q[1];
rz(2.887602) q[1];
x q[2];
rz(-2.0207094) q[3];
sx q[3];
rz(-0.16928798) q[3];
sx q[3];
rz(-1.560488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70904237) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(0.77152983) q[2];
rz(-0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(-2.5468723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58525697) q[0];
sx q[0];
rz(-0.62194967) q[0];
sx q[0];
rz(2.0167895) q[0];
rz(-1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(-0.76464701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2128396) q[0];
sx q[0];
rz(-3.0055586) q[0];
sx q[0];
rz(0.68599756) q[0];
rz(-pi) q[1];
rz(-2.6420399) q[2];
sx q[2];
rz(-0.95328125) q[2];
sx q[2];
rz(-2.2411186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5950039) q[1];
sx q[1];
rz(-1.6975132) q[1];
sx q[1];
rz(-1.3664043) q[1];
rz(-1.9493136) q[3];
sx q[3];
rz(-2.0928185) q[3];
sx q[3];
rz(1.4498364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(-2.1616705) q[3];
sx q[3];
rz(-2.0082974) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205263) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(3.0038515) q[2];
sx q[2];
rz(-1.6888035) q[2];
sx q[2];
rz(2.9464108) q[2];
rz(-1.408314) q[3];
sx q[3];
rz(-1.8058895) q[3];
sx q[3];
rz(-1.4598082) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
