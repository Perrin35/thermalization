OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(-1.4821948) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1971561) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(-0.42874254) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(-0.53127015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2130148) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(1.0282474) q[1];
x q[2];
rz(-0.10098884) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(-1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-0.026219333) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.11461) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(-1.5747889) q[0];
rz(-2.1396779) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(-0.98904726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.023167921) q[1];
sx q[1];
rz(-1.8903362) q[1];
sx q[1];
rz(2.097514) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24641896) q[3];
sx q[3];
rz(-1.2296457) q[3];
sx q[3];
rz(2.6951172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(0.13452402) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(0.55999666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0021734) q[0];
sx q[0];
rz(-1.9529441) q[0];
sx q[0];
rz(2.9234773) q[0];
rz(-pi) q[1];
rz(-1.5921002) q[2];
sx q[2];
rz(-1.2676123) q[2];
sx q[2];
rz(1.6456749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9251688) q[1];
sx q[1];
rz(-2.7270626) q[1];
sx q[1];
rz(0.51149909) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47645724) q[3];
sx q[3];
rz(-2.0072848) q[3];
sx q[3];
rz(-2.1863294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(2.8570783) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38628681) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(-2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5721333) q[2];
sx q[2];
rz(-0.7664116) q[2];
sx q[2];
rz(0.02174755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5609834) q[1];
sx q[1];
rz(-1.0099851) q[1];
sx q[1];
rz(2.9210655) q[1];
x q[2];
rz(-1.8989765) q[3];
sx q[3];
rz(-1.4644074) q[3];
sx q[3];
rz(2.1387517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(-0.38468012) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.4543021) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41243991) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(-0.57106437) q[0];
rz(-0.16634059) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(1.2321842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2424803) q[1];
sx q[1];
rz(-2.50933) q[1];
sx q[1];
rz(2.2805023) q[1];
rz(0.71367587) q[3];
sx q[3];
rz(-1.5034961) q[3];
sx q[3];
rz(-2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(2.5352535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0814708) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(1.5541374) q[0];
x q[1];
rz(-2.7251284) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(3.0481899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50748435) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(-0.21168153) q[1];
rz(-pi) q[2];
rz(-0.074999853) q[3];
sx q[3];
rz(-1.3538176) q[3];
sx q[3];
rz(-1.0523588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63885826) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(1.139337) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.5320324) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-0.79024822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797487) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(2.5248812) q[0];
x q[1];
rz(0.054025606) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(0.33188785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(-0.40562628) q[1];
rz(-1.1931476) q[3];
sx q[3];
rz(-0.95818633) q[3];
sx q[3];
rz(0.8876422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097276) q[0];
sx q[0];
rz(-2.1301529) q[0];
sx q[0];
rz(3.1064242) q[0];
rz(-2.6128204) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(0.96226529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(0.1235048) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96447585) q[3];
sx q[3];
rz(-1.3250321) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(-0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(2.5833599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24960625) q[0];
sx q[0];
rz(-0.59392801) q[0];
sx q[0];
rz(2.3114165) q[0];
x q[1];
rz(-0.054106601) q[2];
sx q[2];
rz(-1.4634842) q[2];
sx q[2];
rz(0.12003128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25591125) q[1];
sx q[1];
rz(-0.7464039) q[1];
sx q[1];
rz(1.6649151) q[1];
rz(-pi) q[2];
rz(-1.3853119) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(-0.51717796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.5378392) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1062746) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(1.6465228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043838219) q[2];
sx q[2];
rz(-0.99928108) q[2];
sx q[2];
rz(2.8457355) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.938898) q[1];
sx q[1];
rz(-1.885186) q[1];
sx q[1];
rz(2.7640192) q[1];
rz(-pi) q[2];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(0.45104879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(2.8722897) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.6745463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(3.1303828) q[2];
sx q[2];
rz(-1.3080454) q[2];
sx q[2];
rz(-1.0946454) q[2];
rz(1.8105103) q[3];
sx q[3];
rz(-2.3532383) q[3];
sx q[3];
rz(-0.91010703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
