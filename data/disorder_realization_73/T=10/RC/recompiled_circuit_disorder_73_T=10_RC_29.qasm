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
rz(2.4174262) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1278348) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(-0.6291997) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6720812) q[2];
sx q[2];
rz(-1.9874319) q[2];
sx q[2];
rz(0.11726221) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29014978) q[1];
sx q[1];
rz(-2.5622497) q[1];
sx q[1];
rz(1.9860553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4943487) q[3];
sx q[3];
rz(-0.41562286) q[3];
sx q[3];
rz(1.7474183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.589754) q[2];
sx q[2];
rz(-1.4171615) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-0.63175732) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(-1.7858645) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512017) q[0];
sx q[0];
rz(-2.8848007) q[0];
sx q[0];
rz(-1.4635565) q[0];
rz(-pi) q[1];
rz(1.9594876) q[2];
sx q[2];
rz(-1.9027862) q[2];
sx q[2];
rz(1.7322844) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4888549) q[1];
sx q[1];
rz(-1.6994209) q[1];
sx q[1];
rz(-1.1210404) q[1];
rz(-pi) q[2];
rz(-1.5242819) q[3];
sx q[3];
rz(-2.1624613) q[3];
sx q[3];
rz(3.068963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(-2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8963985) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(-1.3348745) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.867884) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-2.6991978) q[0];
sx q[0];
rz(0.94985234) q[0];
x q[1];
rz(-1.4726228) q[2];
sx q[2];
rz(-1.86491) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0395567) q[1];
sx q[1];
rz(-1.9196871) q[1];
sx q[1];
rz(-2.4400913) q[1];
rz(-pi) q[2];
rz(-0.55176118) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(1.7337652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(2.188142) q[2];
rz(0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-0.24818534) q[0];
rz(1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1124681) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(2.8549457) q[0];
rz(-0.3481625) q[2];
sx q[2];
rz(-2.2194038) q[2];
sx q[2];
rz(-1.6096889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3809507) q[1];
sx q[1];
rz(-1.6886854) q[1];
sx q[1];
rz(-1.2537434) q[1];
rz(-0.62319237) q[3];
sx q[3];
rz(-2.8487301) q[3];
sx q[3];
rz(-2.5588536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(1.3624396) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.978925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8602596) q[0];
sx q[0];
rz(-1.6006002) q[0];
sx q[0];
rz(1.4593065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9722399) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(-0.40118518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77046493) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(-2.0626555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.053247) q[3];
sx q[3];
rz(-1.9540817) q[3];
sx q[3];
rz(-0.41107197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(2.999372) q[2];
rz(-0.90406117) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-2.8335559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7361569) q[0];
sx q[0];
rz(-1.476256) q[0];
sx q[0];
rz(-0.12673881) q[0];
rz(0.61200895) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(-2.3305364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(-1.874079) q[1];
rz(-pi) q[2];
rz(-2.8940053) q[3];
sx q[3];
rz(-2.7831804) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(-1.1479088) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797453) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(-3.1107483) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(2.470509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4422465) q[0];
sx q[0];
rz(-1.6244349) q[0];
sx q[0];
rz(-1.5357114) q[0];
x q[1];
rz(1.7332156) q[2];
sx q[2];
rz(-0.186609) q[2];
sx q[2];
rz(-0.23255541) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1783501) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(-0.11022719) q[1];
rz(-0.047366553) q[3];
sx q[3];
rz(-0.8960552) q[3];
sx q[3];
rz(-0.55344068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-2.7291765) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0359356) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(-2.1980594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2799759) q[2];
sx q[2];
rz(-1.78252) q[2];
sx q[2];
rz(2.784563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7653212) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(-0.18860753) q[1];
rz(1.8230209) q[3];
sx q[3];
rz(-0.54402292) q[3];
sx q[3];
rz(-0.23214425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2010487) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(-1.3930901) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.4896726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.888436) q[0];
sx q[0];
rz(-0.90421593) q[0];
sx q[0];
rz(-2.1783834) q[0];
x q[1];
rz(1.4719047) q[2];
sx q[2];
rz(-0.47404587) q[2];
sx q[2];
rz(-0.74741077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1676294) q[1];
sx q[1];
rz(-1.3988252) q[1];
sx q[1];
rz(-0.19584882) q[1];
rz(-pi) q[2];
rz(0.49976607) q[3];
sx q[3];
rz(-1.0904113) q[3];
sx q[3];
rz(1.2479316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5513409) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518236) q[0];
sx q[0];
rz(-0.71634403) q[0];
sx q[0];
rz(0.96869529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5006127) q[2];
sx q[2];
rz(-1.3314221) q[2];
sx q[2];
rz(-2.0812876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30565572) q[1];
sx q[1];
rz(-0.61969295) q[1];
sx q[1];
rz(-1.6395007) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57659984) q[3];
sx q[3];
rz(-1.7017662) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(3.0030737) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(1.7759454) q[3];
sx q[3];
rz(-0.57153265) q[3];
sx q[3];
rz(2.3415757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
