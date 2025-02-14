OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(-1.1298236) q[0];
sx q[0];
rz(3.1273754) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(-2.3925233) q[1];
sx q[1];
rz(0.60325375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39122836) q[0];
sx q[0];
rz(-0.75699139) q[0];
sx q[0];
rz(-2.6177177) q[0];
rz(-0.91602625) q[2];
sx q[2];
rz(-1.6189304) q[2];
sx q[2];
rz(2.2544238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7426213) q[1];
sx q[1];
rz(-1.7265897) q[1];
sx q[1];
rz(0.69424082) q[1];
rz(2.2774062) q[3];
sx q[3];
rz(-1.4281311) q[3];
sx q[3];
rz(-0.29970009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3732036) q[2];
sx q[2];
rz(-1.5378121) q[2];
sx q[2];
rz(-2.0453889) q[2];
rz(-1.0783892) q[3];
sx q[3];
rz(-2.6069141) q[3];
sx q[3];
rz(-0.53511846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32644367) q[0];
sx q[0];
rz(-1.0468227) q[0];
sx q[0];
rz(0.4775508) q[0];
rz(2.1304456) q[1];
sx q[1];
rz(-0.54713455) q[1];
sx q[1];
rz(1.0844885) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835533) q[0];
sx q[0];
rz(-0.38738444) q[0];
sx q[0];
rz(-1.256146) q[0];
rz(-pi) q[1];
rz(-1.6376375) q[2];
sx q[2];
rz(-2.7324841) q[2];
sx q[2];
rz(2.6928549) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66531679) q[1];
sx q[1];
rz(-1.3937104) q[1];
sx q[1];
rz(-0.55426532) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7803187) q[3];
sx q[3];
rz(-0.84224904) q[3];
sx q[3];
rz(-1.1948578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61887211) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(-2.0514533) q[2];
rz(2.5777396) q[3];
sx q[3];
rz(-2.4520051) q[3];
sx q[3];
rz(3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801333) q[0];
sx q[0];
rz(-3.1205803) q[0];
sx q[0];
rz(-1.6047961) q[0];
rz(1.3179294) q[1];
sx q[1];
rz(-2.047796) q[1];
sx q[1];
rz(0.39753786) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73765342) q[0];
sx q[0];
rz(-1.7884457) q[0];
sx q[0];
rz(-2.2277839) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8454281) q[2];
sx q[2];
rz(-2.1284747) q[2];
sx q[2];
rz(-2.9147343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6315115) q[1];
sx q[1];
rz(-0.95247522) q[1];
sx q[1];
rz(0.31429283) q[1];
rz(-pi) q[2];
rz(-1.0646754) q[3];
sx q[3];
rz(-2.0398413) q[3];
sx q[3];
rz(-2.7739547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2354108) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(-2.5664491) q[2];
rz(-2.4681674) q[3];
sx q[3];
rz(-1.5410825) q[3];
sx q[3];
rz(-1.5448236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49066082) q[0];
sx q[0];
rz(-2.0199825) q[0];
sx q[0];
rz(-2.2326873) q[0];
rz(0.017223651) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-3.1294894) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229268) q[0];
sx q[0];
rz(-1.2483435) q[0];
sx q[0];
rz(-0.61005436) q[0];
x q[1];
rz(2.4922089) q[2];
sx q[2];
rz(-2.0133791) q[2];
sx q[2];
rz(-1.8049406) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.035701) q[1];
sx q[1];
rz(-1.3451223) q[1];
sx q[1];
rz(-2.1848328) q[1];
rz(-0.34194591) q[3];
sx q[3];
rz(-0.13445915) q[3];
sx q[3];
rz(-2.16914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0688613) q[2];
sx q[2];
rz(-0.27421811) q[2];
sx q[2];
rz(0.38632986) q[2];
rz(-2.7522411) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(-1.0079591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35578457) q[0];
sx q[0];
rz(-0.72440994) q[0];
sx q[0];
rz(2.8177148) q[0];
rz(-2.7549699) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(1.5547543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6232658) q[0];
sx q[0];
rz(-0.73438209) q[0];
sx q[0];
rz(3.0452888) q[0];
rz(-pi) q[1];
rz(0.37347157) q[2];
sx q[2];
rz(-1.0658776) q[2];
sx q[2];
rz(-0.89044705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2302479) q[1];
sx q[1];
rz(-0.94083386) q[1];
sx q[1];
rz(2.1542284) q[1];
rz(-pi) q[2];
rz(-2.6324248) q[3];
sx q[3];
rz(-0.8664136) q[3];
sx q[3];
rz(2.1455163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(-3.0774934) q[2];
rz(0.60424232) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(1.2324246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9391249) q[0];
sx q[0];
rz(-2.5305643) q[0];
sx q[0];
rz(2.2623999) q[0];
rz(0.35036707) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(2.9894357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50987303) q[0];
sx q[0];
rz(-1.3001967) q[0];
sx q[0];
rz(0.71978777) q[0];
rz(-pi) q[1];
rz(-2.9486676) q[2];
sx q[2];
rz(-2.0461651) q[2];
sx q[2];
rz(-2.1383291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81318839) q[1];
sx q[1];
rz(-1.1394573) q[1];
sx q[1];
rz(-1.3479517) q[1];
rz(-pi) q[2];
x q[2];
rz(2.426126) q[3];
sx q[3];
rz(-2.5757534) q[3];
sx q[3];
rz(0.51763207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0396314) q[2];
sx q[2];
rz(-0.32932082) q[2];
sx q[2];
rz(2.9015818) q[2];
rz(0.65252423) q[3];
sx q[3];
rz(-2.0034761) q[3];
sx q[3];
rz(1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(2.2834593) q[0];
rz(-1.3872604) q[1];
sx q[1];
rz(-2.6871082) q[1];
sx q[1];
rz(2.8087356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3368206) q[0];
sx q[0];
rz(-0.66786486) q[0];
sx q[0];
rz(2.0940381) q[0];
rz(-pi) q[1];
rz(0.24999745) q[2];
sx q[2];
rz(-1.6277024) q[2];
sx q[2];
rz(3.0661063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0197647) q[1];
sx q[1];
rz(-0.47224829) q[1];
sx q[1];
rz(2.6761503) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4868785) q[3];
sx q[3];
rz(-1.2867478) q[3];
sx q[3];
rz(-2.9971307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7994999) q[2];
sx q[2];
rz(-2.9436593) q[2];
sx q[2];
rz(-1.0510772) q[2];
rz(-1.6445232) q[3];
sx q[3];
rz(-1.9688508) q[3];
sx q[3];
rz(-0.79819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77338162) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(-2.2494702) q[0];
rz(-0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(-0.77286744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5839856) q[0];
sx q[0];
rz(-2.0933042) q[0];
sx q[0];
rz(0.64675348) q[0];
rz(-pi) q[1];
rz(2.4382227) q[2];
sx q[2];
rz(-1.7711319) q[2];
sx q[2];
rz(0.021066152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9443431) q[1];
sx q[1];
rz(-2.5641003) q[1];
sx q[1];
rz(0.61402278) q[1];
rz(-pi) q[2];
rz(2.3691133) q[3];
sx q[3];
rz(-2.8101343) q[3];
sx q[3];
rz(-0.14642388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.2246776) q[2];
sx q[2];
rz(2.4490855) q[2];
rz(0.45037371) q[3];
sx q[3];
rz(-1.2053442) q[3];
sx q[3];
rz(1.5041806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5699128) q[0];
sx q[0];
rz(-2.3501861) q[0];
sx q[0];
rz(-2.1929542) q[0];
rz(1.0665077) q[1];
sx q[1];
rz(-2.152161) q[1];
sx q[1];
rz(1.271064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1656837) q[0];
sx q[0];
rz(-1.6635487) q[0];
sx q[0];
rz(1.776987) q[0];
rz(-pi) q[1];
rz(-1.6048715) q[2];
sx q[2];
rz(-1.5487031) q[2];
sx q[2];
rz(1.3349229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3349325) q[1];
sx q[1];
rz(-0.94881034) q[1];
sx q[1];
rz(-2.9750664) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27008121) q[3];
sx q[3];
rz(-0.45680667) q[3];
sx q[3];
rz(2.1169259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54242414) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(0.14031169) q[2];
rz(-3.1091651) q[3];
sx q[3];
rz(-0.92493886) q[3];
sx q[3];
rz(1.2062581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5133544) q[0];
sx q[0];
rz(-1.6980549) q[0];
sx q[0];
rz(2.3640609) q[0];
rz(-2.2041722) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(-2.6984528) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86771482) q[0];
sx q[0];
rz(-1.4860503) q[0];
sx q[0];
rz(1.4895208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1681917) q[2];
sx q[2];
rz(-1.294916) q[2];
sx q[2];
rz(-1.2593002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3819246) q[1];
sx q[1];
rz(-0.72205359) q[1];
sx q[1];
rz(-0.51314919) q[1];
rz(-pi) q[2];
rz(-0.16297131) q[3];
sx q[3];
rz(-1.2633557) q[3];
sx q[3];
rz(-0.93443894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-3.013986) q[2];
rz(-2.8429032) q[3];
sx q[3];
rz(-1.1726215) q[3];
sx q[3];
rz(2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61113197) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(-0.61027377) q[1];
sx q[1];
rz(-0.88519575) q[1];
sx q[1];
rz(-2.0559678) q[1];
rz(0.064432003) q[2];
sx q[2];
rz(-2.2885135) q[2];
sx q[2];
rz(2.0533333) q[2];
rz(0.21214938) q[3];
sx q[3];
rz(-1.129088) q[3];
sx q[3];
rz(-1.3442985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
