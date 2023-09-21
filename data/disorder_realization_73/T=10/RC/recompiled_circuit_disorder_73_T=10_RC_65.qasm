OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(-0.72416645) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35207507) q[0];
sx q[0];
rz(-1.8624458) q[0];
sx q[0];
rz(-1.7910936) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22462331) q[2];
sx q[2];
rz(-0.42806872) q[2];
sx q[2];
rz(-3.012804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6341056) q[1];
sx q[1];
rz(-1.7934985) q[1];
sx q[1];
rz(2.1102064) q[1];
rz(-pi) q[2];
rz(-1.9853398) q[3];
sx q[3];
rz(-1.539955) q[3];
sx q[3];
rz(0.24658345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(-0.12456482) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403554) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(-3.113494) q[0];
rz(-1.9594876) q[2];
sx q[2];
rz(-1.9027862) q[2];
sx q[2];
rz(1.4093083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6527378) q[1];
sx q[1];
rz(-1.6994209) q[1];
sx q[1];
rz(1.1210404) q[1];
rz(-pi) q[2];
rz(3.0725067) q[3];
sx q[3];
rz(-2.54832) q[3];
sx q[3];
rz(0.010635076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0624861) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(0.50659531) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(-0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(2.2431592) q[0];
rz(-1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.2737087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(-0.94985234) q[0];
x q[1];
rz(2.8461371) q[2];
sx q[2];
rz(-1.4768512) q[2];
sx q[2];
rz(0.23677793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18759218) q[1];
sx q[1];
rz(-2.2224269) q[1];
sx q[1];
rz(1.1263532) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55176118) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0597824) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(-2.10363) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(-3.0674556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.603133) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(2.0752226) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2494227) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-0.25472578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84872765) q[1];
sx q[1];
rz(-1.25602) q[1];
sx q[1];
rz(-3.0175812) q[1];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(-3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8903824) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(1.779153) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(1.978925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8602596) q[0];
sx q[0];
rz(-1.6006002) q[0];
sx q[0];
rz(1.6822862) q[0];
rz(-pi) q[1];
rz(-1.8936929) q[2];
sx q[2];
rz(-1.731551) q[2];
sx q[2];
rz(-1.9184743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77046493) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(1.0789372) q[1];
x q[2];
rz(2.707162) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.7720951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(-0.18946762) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(1.4676771) q[0];
rz(0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(0.30803672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98826) q[0];
sx q[0];
rz(-1.4446265) q[0];
sx q[0];
rz(-1.6660965) q[0];
x q[1];
rz(2.1790444) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(0.42524291) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.249237) q[1];
sx q[1];
rz(-0.71214572) q[1];
sx q[1];
rz(-1.2675136) q[1];
rz(-2.7932348) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-0.1299783) q[0];
rz(-3.1107483) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-2.470509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86258793) q[0];
sx q[0];
rz(-0.064084856) q[0];
sx q[0];
rz(2.5628753) q[0];
x q[1];
rz(-1.3865878) q[2];
sx q[2];
rz(-1.5407908) q[2];
sx q[2];
rz(-1.9629994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.259686) q[1];
sx q[1];
rz(-2.7555008) q[1];
sx q[1];
rz(1.2950456) q[1];
rz(-pi) q[2];
rz(-1.6299134) q[3];
sx q[3];
rz(-0.67614188) q[3];
sx q[3];
rz(-2.6638871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-2.7291765) q[0];
rz(1.4498129) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(1.9746045) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0359356) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(0.94353326) q[0];
x q[1];
rz(1.251986) q[2];
sx q[2];
rz(-0.73482162) q[2];
sx q[2];
rz(-2.1679945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3762714) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(-2.9529851) q[1];
x q[2];
rz(-2.1006881) q[3];
sx q[3];
rz(-1.7003254) q[3];
sx q[3];
rz(-2.0199752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-0.92700672) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.6519201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.888436) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(2.1783834) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4719047) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(2.3941819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1676294) q[1];
sx q[1];
rz(-1.7427674) q[1];
sx q[1];
rz(-0.19584882) q[1];
x q[2];
rz(-2.6418266) q[3];
sx q[3];
rz(-1.0904113) q[3];
sx q[3];
rz(1.2479316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(-1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9592181) q[0];
sx q[0];
rz(-1.1897414) q[0];
sx q[0];
rz(-2.1931838) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9016501) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(-0.49382526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22132561) q[1];
sx q[1];
rz(-2.1888071) q[1];
sx q[1];
rz(3.0926535) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5649928) q[3];
sx q[3];
rz(-1.7017662) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(-0.74404136) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-2.3241282) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(1.6384009) q[2];
sx q[2];
rz(-1.1195782) q[2];
sx q[2];
rz(-1.3200214) q[2];
rz(-1.0088624) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
