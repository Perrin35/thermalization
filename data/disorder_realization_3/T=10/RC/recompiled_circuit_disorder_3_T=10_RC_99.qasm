OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1026693) q[0];
sx q[0];
rz(-1.9063213) q[0];
sx q[0];
rz(-1.194792) q[0];
rz(-pi) q[1];
rz(-2.1404999) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(-0.18730883) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2990883) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(-1.7249785) q[1];
rz(-pi) q[2];
rz(-1.4685417) q[3];
sx q[3];
rz(-1.8342606) q[3];
sx q[3];
rz(-1.2276358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182506) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(-2.0204861) q[0];
rz(-pi) q[1];
rz(-1.1621446) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(-2.0603927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0961699) q[1];
sx q[1];
rz(-1.759937) q[1];
sx q[1];
rz(-2.0778836) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5856421) q[3];
sx q[3];
rz(-0.9405989) q[3];
sx q[3];
rz(0.48660183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(-2.4754751) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(2.8853436) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51107823) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(-1.2468673) q[0];
x q[1];
rz(-1.0923549) q[2];
sx q[2];
rz(-1.97176) q[2];
sx q[2];
rz(-1.5329597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5445404) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(0.87045963) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17685299) q[3];
sx q[3];
rz(-1.6079418) q[3];
sx q[3];
rz(1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-0.57470542) q[2];
rz(-1.3556708) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.2333966) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(-0.5831232) q[0];
x q[1];
rz(3.0078997) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(2.0267817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(2.494032) q[1];
x q[2];
rz(1.6963523) q[3];
sx q[3];
rz(-1.7896277) q[3];
sx q[3];
rz(-2.9343176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3691468) q[0];
sx q[0];
rz(-1.4506842) q[0];
sx q[0];
rz(-0.0025047501) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30349489) q[2];
sx q[2];
rz(-1.3611088) q[2];
sx q[2];
rz(1.4022624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.55802) q[1];
sx q[1];
rz(-2.1961374) q[1];
sx q[1];
rz(3.0261092) q[1];
x q[2];
rz(-0.4425211) q[3];
sx q[3];
rz(-1.2708775) q[3];
sx q[3];
rz(0.72344852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698343) q[0];
sx q[0];
rz(-1.7239128) q[0];
sx q[0];
rz(-1.7437177) q[0];
x q[1];
rz(2.3577865) q[2];
sx q[2];
rz(-1.1597826) q[2];
sx q[2];
rz(-2.59936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97092123) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(1.3794273) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.4962713) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(0.62136674) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646274) q[0];
sx q[0];
rz(-1.0210751) q[0];
sx q[0];
rz(-1.2697898) q[0];
rz(0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-2.0932587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.36564) q[1];
sx q[1];
rz(-1.0675634) q[1];
sx q[1];
rz(-0.788049) q[1];
rz(1.9637397) q[3];
sx q[3];
rz(-1.8624536) q[3];
sx q[3];
rz(1.1100811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(-0.41637862) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(2.7767048) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.5023124) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84425981) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(0.83606007) q[2];
sx q[2];
rz(-1.2088472) q[2];
sx q[2];
rz(-1.8723633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2316206) q[1];
sx q[1];
rz(-0.93712229) q[1];
sx q[1];
rz(2.2909067) q[1];
rz(-pi) q[2];
rz(0.092076093) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(0.74842194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(0.41697821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96796658) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(-3.0941512) q[0];
x q[1];
rz(2.142799) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(-1.4243766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7855362) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(-2.9836125) q[1];
rz(-pi) q[2];
rz(-3.0005089) q[3];
sx q[3];
rz(-1.4420296) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.8189925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(1.9508719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(2.3479455) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(-0.91099693) q[1];
x q[2];
rz(-2.2262276) q[3];
sx q[3];
rz(-1.6928821) q[3];
sx q[3];
rz(-1.9742427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.017756391) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(-0.090311269) q[3];
sx q[3];
rz(-1.8009381) q[3];
sx q[3];
rz(1.2931852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];