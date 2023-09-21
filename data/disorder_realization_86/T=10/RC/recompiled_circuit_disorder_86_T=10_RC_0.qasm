OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(-1.7295184) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(-1.4821948) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5858551) q[0];
sx q[0];
rz(-1.9940388) q[0];
sx q[0];
rz(-1.4002355) q[0];
rz(-pi) q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-2.854752) q[2];
sx q[2];
rz(-1.6490205) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(-0.4835101) q[1];
rz(-1.0079908) q[3];
sx q[3];
rz(-1.6562914) q[3];
sx q[3];
rz(2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(-0.95430294) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-0.96347934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.11461) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(-1.5668037) q[0];
x q[1];
rz(2.1396779) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(0.98904726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7746437) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(-0.36555396) q[1];
rz(1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-3.0070686) q[2];
rz(-2.3965805) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(2.581596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0021734) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(2.9234773) q[0];
rz(-pi) q[1];
rz(-0.30324869) q[2];
sx q[2];
rz(-1.5504642) q[2];
sx q[2];
rz(3.0603527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9251688) q[1];
sx q[1];
rz(-0.41453002) q[1];
sx q[1];
rz(2.6300936) q[1];
rz(-pi) q[2];
rz(0.47645724) q[3];
sx q[3];
rz(-2.0072848) q[3];
sx q[3];
rz(-2.1863294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(-2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383706) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(-0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38628681) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(1.132071) q[0];
x q[1];
rz(-0.80438517) q[2];
sx q[2];
rz(-1.569869) q[2];
sx q[2];
rz(-1.550012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.270077) q[1];
sx q[1];
rz(-1.3844826) q[1];
sx q[1];
rz(-0.99888505) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11233791) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(1.4543021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(-2.9105913) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-0.2968266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2192229) q[0];
sx q[0];
rz(-2.5693359) q[0];
sx q[0];
rz(3.0692283) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7237687) q[2];
sx q[2];
rz(-2.3077871) q[2];
sx q[2];
rz(-2.134915) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2424803) q[1];
sx q[1];
rz(-2.50933) q[1];
sx q[1];
rz(-2.2805023) q[1];
rz(-pi) q[2];
rz(3.0389901) q[3];
sx q[3];
rz(-2.4253064) q[3];
sx q[3];
rz(-0.41075452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(-1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
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
rz(-1.5874552) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41646429) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(3.0481899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9669173) q[1];
sx q[1];
rz(-1.7511909) q[1];
sx q[1];
rz(-1.0134407) q[1];
x q[2];
rz(-1.7883676) q[3];
sx q[3];
rz(-1.6440344) q[3];
sx q[3];
rz(0.53461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-1.139337) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(2.3513444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7471874) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(1.2140843) q[0];
rz(0.054025606) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(0.33188785) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0223169) q[1];
sx q[1];
rz(-1.1806618) q[1];
sx q[1];
rz(-1.2783865) q[1];
x q[2];
rz(2.6584133) q[3];
sx q[3];
rz(-0.70671591) q[3];
sx q[3];
rz(0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-1.4247894) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(-2.5324902) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(-1.3887127) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213288) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(1.0111615) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52877229) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(-0.96226529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70506239) q[1];
sx q[1];
rz(-3.0490746) q[1];
sx q[1];
rz(-0.1235048) q[1];
rz(-0.29626131) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(-2.0010009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9528815) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(-0.88225538) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6728014) q[0];
sx q[0];
rz(-1.9580012) q[0];
sx q[0];
rz(1.1084491) q[0];
rz(-pi) q[1];
rz(1.105537) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(-0.82660316) q[1];
rz(-pi) q[2];
rz(-1.7562808) q[3];
sx q[3];
rz(-2.0501325) q[3];
sx q[3];
rz(2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-2.9294087) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50487173) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16738811) q[0];
sx q[0];
rz(-2.5295969) q[0];
sx q[0];
rz(-0.1083072) q[0];
rz(-pi) q[1];
rz(0.9988437) q[2];
sx q[2];
rz(-1.5339282) q[2];
sx q[2];
rz(-1.2986623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0303505) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-0.72279795) q[1];
rz(-pi) q[2];
rz(-0.75533112) q[3];
sx q[3];
rz(-0.33077251) q[3];
sx q[3];
rz(1.2942693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-2.0555029) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257278) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.4670463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(1.3080296) q[2];
sx q[2];
rz(-1.5816214) q[2];
sx q[2];
rz(-2.66253) q[2];
rz(-2.3446463) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
