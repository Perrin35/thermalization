OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46086754) q[0];
sx q[0];
rz(-1.4493677) q[0];
sx q[0];
rz(0.073343883) q[0];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(-1.4717799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.764896) q[1];
sx q[1];
rz(-2.2590859) q[1];
sx q[1];
rz(0.98801686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7226726) q[3];
sx q[3];
rz(-1.6994611) q[3];
sx q[3];
rz(-3.0959689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(2.4689891) q[2];
rz(2.9721695) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(0.16054343) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-0.28796089) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32591336) q[0];
sx q[0];
rz(-0.26198146) q[0];
sx q[0];
rz(2.322305) q[0];
rz(-2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(0.39272768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8781232) q[1];
sx q[1];
rz(-2.1681004) q[1];
sx q[1];
rz(2.0720481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5324627) q[3];
sx q[3];
rz(-1.9209071) q[3];
sx q[3];
rz(2.4083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8163452) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(2.8741799) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(0.95169383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0824453) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(-3.0491769) q[0];
rz(2.0864262) q[2];
sx q[2];
rz(-1.3736758) q[2];
sx q[2];
rz(2.4899763) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0905076) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(-2.3329263) q[1];
x q[2];
rz(0.53593105) q[3];
sx q[3];
rz(-1.0169573) q[3];
sx q[3];
rz(1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(-2.7974131) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-0.27799806) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(1.8970867) q[0];
rz(-0.1291153) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(-0.37277645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059797212) q[0];
sx q[0];
rz(-0.76515388) q[0];
sx q[0];
rz(-2.9761936) q[0];
rz(-2.3635025) q[2];
sx q[2];
rz(-1.2568682) q[2];
sx q[2];
rz(0.84184605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72153889) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(-2.0481471) q[1];
rz(-pi) q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.8653449) q[3];
sx q[3];
rz(1.095872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-2.2367031) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(1.6632535) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(0.53422654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.8792218) q[0];
sx q[0];
rz(-2.8240859) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0787557) q[2];
sx q[2];
rz(-1.2675708) q[2];
sx q[2];
rz(-2.6360896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6492669) q[1];
sx q[1];
rz(-1.5070792) q[1];
sx q[1];
rz(1.6462412) q[1];
rz(-pi) q[2];
rz(-1.7027431) q[3];
sx q[3];
rz(-0.85634106) q[3];
sx q[3];
rz(0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(-1.2409695) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.7664849) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2415337) q[0];
sx q[0];
rz(-2.7068479) q[0];
sx q[0];
rz(-2.4777806) q[0];
rz(2.3926922) q[2];
sx q[2];
rz(-1.5178711) q[2];
sx q[2];
rz(0.63268328) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75406972) q[1];
sx q[1];
rz(-1.9807528) q[1];
sx q[1];
rz(-1.5974664) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70011219) q[3];
sx q[3];
rz(-0.60986076) q[3];
sx q[3];
rz(-1.3175347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(-0.14870816) q[2];
rz(3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-0.08392863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5556363) q[0];
sx q[0];
rz(-1.7935828) q[0];
sx q[0];
rz(-1.0046093) q[0];
rz(-pi) q[1];
rz(-1.3533808) q[2];
sx q[2];
rz(-1.2251717) q[2];
sx q[2];
rz(-2.2923922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2578656) q[1];
sx q[1];
rz(-0.95974937) q[1];
sx q[1];
rz(-2.3962767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84079068) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(-1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-2.596358) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(0.30495131) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298744) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(-0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0154008) q[0];
sx q[0];
rz(-1.3500824) q[0];
sx q[0];
rz(1.4908233) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85407599) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(2.9727109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7860124) q[1];
sx q[1];
rz(-1.542172) q[1];
sx q[1];
rz(-0.37866681) q[1];
rz(-pi) q[2];
rz(-1.5787015) q[3];
sx q[3];
rz(-2.2369011) q[3];
sx q[3];
rz(2.283309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(0.61532414) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(-2.5659134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262779) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(-2.4995575) q[0];
rz(-pi) q[1];
rz(-1.3747146) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(-0.17734222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5902192) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(1.3361841) q[1];
x q[2];
rz(1.2647763) q[3];
sx q[3];
rz(-2.5796606) q[3];
sx q[3];
rz(2.7276873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(0.021961948) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-2.172487) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983738) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(1.6427342) q[0];
rz(-1*pi/15) q[2];
sx q[2];
rz(-2.3839715) q[2];
sx q[2];
rz(0.19291887) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5574054) q[1];
sx q[1];
rz(-2.4483042) q[1];
sx q[1];
rz(2.8597945) q[1];
x q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(-2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071028) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(-0.044152505) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(2.2864441) q[2];
sx q[2];
rz(-0.48968857) q[2];
sx q[2];
rz(-0.80077632) q[2];
rz(0.45331656) q[3];
sx q[3];
rz(-1.286187) q[3];
sx q[3];
rz(-1.5766531) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
