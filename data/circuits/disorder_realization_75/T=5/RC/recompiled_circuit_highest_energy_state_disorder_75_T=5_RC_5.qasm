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
rz(1.3156112) q[0];
sx q[0];
rz(5.4551107) q[0];
sx q[0];
rz(8.5760737) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(0.42660776) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72683452) q[0];
sx q[0];
rz(-1.6462741) q[0];
sx q[0];
rz(0.10087905) q[0];
rz(1.4584862) q[2];
sx q[2];
rz(-1.6717615) q[2];
sx q[2];
rz(1.2580521) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5004896) q[1];
sx q[1];
rz(-1.6436542) q[1];
sx q[1];
rz(-0.44568731) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1900666) q[3];
sx q[3];
rz(-2.2689156) q[3];
sx q[3];
rz(-1.2941456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9244869) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(-3.1021049) q[2];
rz(2.110179) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(1.903681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6295488) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(1.398983) q[0];
rz(-1.6276739) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(-1.4167851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607781) q[0];
sx q[0];
rz(-2.3153458) q[0];
sx q[0];
rz(-1.2673402) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.879989) q[2];
sx q[2];
rz(-0.79643476) q[2];
sx q[2];
rz(-1.6557882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7287711) q[1];
sx q[1];
rz(-1.1755474) q[1];
sx q[1];
rz(1.4229619) q[1];
x q[2];
rz(-1.6187158) q[3];
sx q[3];
rz(-1.1635492) q[3];
sx q[3];
rz(-1.6919235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9713355) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(2.4309168) q[2];
rz(0.014287861) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(-1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9209442) q[0];
sx q[0];
rz(-0.9592239) q[0];
sx q[0];
rz(-0.4162108) q[0];
rz(-1.8572218) q[1];
sx q[1];
rz(-1.6651848) q[1];
sx q[1];
rz(0.31825569) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.672645) q[0];
sx q[0];
rz(-1.8469317) q[0];
sx q[0];
rz(-1.9684845) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2936965) q[2];
sx q[2];
rz(-2.090914) q[2];
sx q[2];
rz(2.0682356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15662312) q[1];
sx q[1];
rz(-0.92069101) q[1];
sx q[1];
rz(1.7414879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8793568) q[3];
sx q[3];
rz(-1.8073544) q[3];
sx q[3];
rz(0.46775888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2661813) q[2];
sx q[2];
rz(-0.77989945) q[2];
sx q[2];
rz(1.8511339) q[2];
rz(-3.087888) q[3];
sx q[3];
rz(-1.3196245) q[3];
sx q[3];
rz(1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5829492) q[0];
sx q[0];
rz(-0.68840331) q[0];
sx q[0];
rz(-1.0097367) q[0];
rz(-2.2263777) q[1];
sx q[1];
rz(-1.6123632) q[1];
sx q[1];
rz(0.42542747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1155419) q[0];
sx q[0];
rz(-2.4049691) q[0];
sx q[0];
rz(1.3213586) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3138626) q[2];
sx q[2];
rz(-1.5166188) q[2];
sx q[2];
rz(2.3119213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76313075) q[1];
sx q[1];
rz(-1.9768856) q[1];
sx q[1];
rz(-0.84880813) q[1];
rz(3.1359659) q[3];
sx q[3];
rz(-1.4079908) q[3];
sx q[3];
rz(0.73782792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2414744) q[2];
sx q[2];
rz(-1.590531) q[2];
sx q[2];
rz(-3.0305064) q[2];
rz(2.9491718) q[3];
sx q[3];
rz(-0.40168732) q[3];
sx q[3];
rz(-2.4470611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7531994) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(-2.2764192) q[0];
rz(0.49258891) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(1.385484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5328952) q[0];
sx q[0];
rz(-2.3430763) q[0];
sx q[0];
rz(0.19143398) q[0];
x q[1];
rz(-2.7215459) q[2];
sx q[2];
rz(-0.85424747) q[2];
sx q[2];
rz(1.6796215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7613449) q[1];
sx q[1];
rz(-1.1780537) q[1];
sx q[1];
rz(2.5775546) q[1];
rz(-pi) q[2];
rz(2.319544) q[3];
sx q[3];
rz(-2.1072142) q[3];
sx q[3];
rz(-2.3182403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8869141) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(-0.706642) q[2];
rz(-2.8142269) q[3];
sx q[3];
rz(-1.8492161) q[3];
sx q[3];
rz(-2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3051598) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(-0.35655546) q[0];
rz(-2.5794079) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(2.1551989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0878746) q[0];
sx q[0];
rz(-2.0656423) q[0];
sx q[0];
rz(-2.3943564) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2012901) q[2];
sx q[2];
rz(-1.8942648) q[2];
sx q[2];
rz(-2.7620014) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2960259) q[1];
sx q[1];
rz(-0.81261501) q[1];
sx q[1];
rz(1.6065099) q[1];
x q[2];
rz(-1.3273193) q[3];
sx q[3];
rz(-1.3946574) q[3];
sx q[3];
rz(1.0276664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76546136) q[2];
sx q[2];
rz(-0.71530801) q[2];
sx q[2];
rz(-0.7005271) q[2];
rz(2.1721407) q[3];
sx q[3];
rz(-1.0337831) q[3];
sx q[3];
rz(-2.2112924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1932061) q[0];
sx q[0];
rz(-0.32873118) q[0];
sx q[0];
rz(2.9902003) q[0];
rz(-1.6311749) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(-1.4600533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23510012) q[0];
sx q[0];
rz(-2.0564046) q[0];
sx q[0];
rz(2.4686345) q[0];
rz(-pi) q[1];
rz(2.1489819) q[2];
sx q[2];
rz(-2.5343203) q[2];
sx q[2];
rz(0.36888514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5826539) q[1];
sx q[1];
rz(-2.451921) q[1];
sx q[1];
rz(0.93780545) q[1];
rz(2.2112956) q[3];
sx q[3];
rz(-0.2411763) q[3];
sx q[3];
rz(-1.1080139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2363362) q[2];
sx q[2];
rz(-0.6826123) q[2];
sx q[2];
rz(-0.35161099) q[2];
rz(0.44175276) q[3];
sx q[3];
rz(-0.92919246) q[3];
sx q[3];
rz(-2.9898804) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0996967) q[0];
sx q[0];
rz(-1.2956887) q[0];
sx q[0];
rz(-1.1610485) q[0];
rz(-0.64758045) q[1];
sx q[1];
rz(-1.7627962) q[1];
sx q[1];
rz(0.82239282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063342) q[0];
sx q[0];
rz(-1.0653235) q[0];
sx q[0];
rz(-2.6837206) q[0];
x q[1];
rz(2.0089125) q[2];
sx q[2];
rz(-0.88399502) q[2];
sx q[2];
rz(-2.7504454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31588512) q[1];
sx q[1];
rz(-1.2824821) q[1];
sx q[1];
rz(-2.350239) q[1];
x q[2];
rz(2.7455198) q[3];
sx q[3];
rz(-1.1669945) q[3];
sx q[3];
rz(1.0084821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2170199) q[2];
sx q[2];
rz(-1.1447039) q[2];
sx q[2];
rz(-1.8355969) q[2];
rz(-2.4235587) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(-3.0927299) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7056535) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(-3.1245681) q[0];
rz(-1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(2.8065525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25680731) q[0];
sx q[0];
rz(-1.5809158) q[0];
sx q[0];
rz(-1.0027893) q[0];
x q[1];
rz(-2.4990988) q[2];
sx q[2];
rz(-1.9064925) q[2];
sx q[2];
rz(3.1106126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6155621) q[1];
sx q[1];
rz(-2.2326734) q[1];
sx q[1];
rz(-0.23961459) q[1];
rz(-pi) q[2];
rz(-2.3565759) q[3];
sx q[3];
rz(-0.64006539) q[3];
sx q[3];
rz(2.9736116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42387858) q[2];
sx q[2];
rz(-2.3953891) q[2];
sx q[2];
rz(0.27900532) q[2];
rz(1.772607) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-2.9947424) q[0];
sx q[0];
rz(-2.3745234) q[0];
rz(2.7686139) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(0.17072089) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6304008) q[0];
sx q[0];
rz(-0.83626473) q[0];
sx q[0];
rz(-1.5319351) q[0];
x q[1];
rz(-2.7266194) q[2];
sx q[2];
rz(-1.7185206) q[2];
sx q[2];
rz(1.8208675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77826277) q[1];
sx q[1];
rz(-2.7365541) q[1];
sx q[1];
rz(-2.2447849) q[1];
rz(-pi) q[2];
rz(1.4108229) q[3];
sx q[3];
rz(-1.108128) q[3];
sx q[3];
rz(1.3079974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79147044) q[2];
sx q[2];
rz(-2.5779724) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(0.9015829) q[3];
sx q[3];
rz(-2.1193347) q[3];
sx q[3];
rz(-0.80488747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881385) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(2.3114655) q[1];
sx q[1];
rz(-1.3273888) q[1];
sx q[1];
rz(-1.8524016) q[1];
rz(2.1813761) q[2];
sx q[2];
rz(-2.001702) q[2];
sx q[2];
rz(-2.3624782) q[2];
rz(-0.64408152) q[3];
sx q[3];
rz(-2.2282176) q[3];
sx q[3];
rz(-0.62716425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
