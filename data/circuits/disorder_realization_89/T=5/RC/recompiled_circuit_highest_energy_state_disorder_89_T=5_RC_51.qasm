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
rz(1.8472449) q[0];
sx q[0];
rz(1.5273153) q[0];
sx q[0];
rz(9.6120678) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(3.6880479) q[1];
sx q[1];
rz(10.772279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2696323) q[0];
sx q[0];
rz(-0.56817164) q[0];
sx q[0];
rz(1.7164036) q[0];
x q[1];
rz(2.7653469) q[2];
sx q[2];
rz(-1.0828138) q[2];
sx q[2];
rz(3.1057292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.627003) q[1];
sx q[1];
rz(-1.5545067) q[1];
sx q[1];
rz(-2.2524627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25661664) q[3];
sx q[3];
rz(-1.3735176) q[3];
sx q[3];
rz(0.2256323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5225141) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(-2.8600233) q[2];
rz(-1.7742026) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(-0.27845964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53681579) q[0];
sx q[0];
rz(-0.29412687) q[0];
sx q[0];
rz(1.7915223) q[0];
rz(0.049909441) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(1.2443589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644783) q[0];
sx q[0];
rz(-1.9196072) q[0];
sx q[0];
rz(-1.6326399) q[0];
rz(-pi) q[1];
rz(0.18644615) q[2];
sx q[2];
rz(-1.0294339) q[2];
sx q[2];
rz(3.0360434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0383953) q[1];
sx q[1];
rz(-1.303325) q[1];
sx q[1];
rz(-2.9385376) q[1];
rz(-pi) q[2];
rz(-0.54173754) q[3];
sx q[3];
rz(-1.1835464) q[3];
sx q[3];
rz(0.70808402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0656978) q[2];
sx q[2];
rz(-1.9588797) q[2];
sx q[2];
rz(0.069570216) q[2];
rz(2.4075497) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(-2.5231584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0760536) q[0];
sx q[0];
rz(-0.19936182) q[0];
sx q[0];
rz(0.25654909) q[0];
rz(-1.8007295) q[1];
sx q[1];
rz(-2.1840265) q[1];
sx q[1];
rz(1.0440913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0437735) q[0];
sx q[0];
rz(-1.7390842) q[0];
sx q[0];
rz(1.2531444) q[0];
x q[1];
rz(1.3568722) q[2];
sx q[2];
rz(-1.7905777) q[2];
sx q[2];
rz(2.8666988) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26790923) q[1];
sx q[1];
rz(-1.2699155) q[1];
sx q[1];
rz(-1.2734702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99756156) q[3];
sx q[3];
rz(-1.4519534) q[3];
sx q[3];
rz(-3.1021743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23036817) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(0.71630859) q[2];
rz(0.45799842) q[3];
sx q[3];
rz(-1.4665946) q[3];
sx q[3];
rz(-0.23076335) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79391795) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(-2.4804261) q[0];
rz(1.2541153) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(1.5987781) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371915) q[0];
sx q[0];
rz(-1.3953623) q[0];
sx q[0];
rz(1.8262922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.132868) q[2];
sx q[2];
rz(-1.9572538) q[2];
sx q[2];
rz(1.5773048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6434272) q[1];
sx q[1];
rz(-1.5076867) q[1];
sx q[1];
rz(-2.4377512) q[1];
rz(-pi) q[2];
rz(-1.0374674) q[3];
sx q[3];
rz(-1.2308443) q[3];
sx q[3];
rz(-1.252591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8864112) q[2];
sx q[2];
rz(-2.7913783) q[2];
sx q[2];
rz(-1.1345471) q[2];
rz(-2.5875731) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968813) q[0];
sx q[0];
rz(-0.0021332707) q[0];
sx q[0];
rz(2.010349) q[0];
rz(3.0834815) q[1];
sx q[1];
rz(-1.2429712) q[1];
sx q[1];
rz(-1.6910472) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3076271) q[0];
sx q[0];
rz(-1.9585607) q[0];
sx q[0];
rz(2.2448756) q[0];
rz(-pi) q[1];
rz(-2.8314231) q[2];
sx q[2];
rz(-0.7855987) q[2];
sx q[2];
rz(0.73332649) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1923545) q[1];
sx q[1];
rz(-2.7325316) q[1];
sx q[1];
rz(2.5109568) q[1];
rz(-pi) q[2];
rz(2.7282894) q[3];
sx q[3];
rz(-2.163475) q[3];
sx q[3];
rz(2.1324219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45945534) q[2];
sx q[2];
rz(-0.61209279) q[2];
sx q[2];
rz(-1.3539782) q[2];
rz(-0.55245429) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(0.49828211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246493) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(1.0981052) q[0];
rz(-1.8662628) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(-1.0414418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324243) q[0];
sx q[0];
rz(-1.0579666) q[0];
sx q[0];
rz(0.70968117) q[0];
rz(-0.11901151) q[2];
sx q[2];
rz(-1.6023984) q[2];
sx q[2];
rz(-1.0234444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2409461) q[1];
sx q[1];
rz(-2.6178872) q[1];
sx q[1];
rz(-1.7745191) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.249635) q[3];
sx q[3];
rz(-0.75157673) q[3];
sx q[3];
rz(-1.0281212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8424524) q[2];
sx q[2];
rz(-1.1536529) q[2];
sx q[2];
rz(3.0843206) q[2];
rz(-2.9412269) q[3];
sx q[3];
rz(-2.5860131) q[3];
sx q[3];
rz(2.949775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054691) q[0];
sx q[0];
rz(-0.64490461) q[0];
sx q[0];
rz(2.7591925) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-2.7303374) q[1];
sx q[1];
rz(1.8866906) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9318274) q[0];
sx q[0];
rz(-1.6269653) q[0];
sx q[0];
rz(2.5302391) q[0];
x q[1];
rz(0.81315984) q[2];
sx q[2];
rz(-0.81612464) q[2];
sx q[2];
rz(-0.94674158) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5211902) q[1];
sx q[1];
rz(-0.56439059) q[1];
sx q[1];
rz(2.317313) q[1];
rz(-0.27886919) q[3];
sx q[3];
rz(-2.3430716) q[3];
sx q[3];
rz(-0.65031933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3249698) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(-0.85186446) q[2];
rz(2.8999117) q[3];
sx q[3];
rz(-1.207374) q[3];
sx q[3];
rz(0.91066256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082751) q[0];
sx q[0];
rz(-2.9589544) q[0];
sx q[0];
rz(2.0727378) q[0];
rz(1.764027) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(0.68881234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4599054) q[0];
sx q[0];
rz(-1.9294159) q[0];
sx q[0];
rz(-1.9214517) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60283349) q[2];
sx q[2];
rz(-1.6844464) q[2];
sx q[2];
rz(1.08504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88676597) q[1];
sx q[1];
rz(-0.51023645) q[1];
sx q[1];
rz(-2.4429823) q[1];
rz(-1.7433002) q[3];
sx q[3];
rz(-2.093708) q[3];
sx q[3];
rz(-0.80181087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0869202) q[2];
sx q[2];
rz(-2.4189147) q[2];
sx q[2];
rz(-1.5642536) q[2];
rz(-2.4857322) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(-2.3724469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5250788) q[0];
sx q[0];
rz(-1.0322796) q[0];
sx q[0];
rz(0.14933625) q[0];
rz(2.0556045) q[1];
sx q[1];
rz(-0.95567742) q[1];
sx q[1];
rz(2.65926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457561) q[0];
sx q[0];
rz(-1.0991503) q[0];
sx q[0];
rz(3.0455941) q[0];
x q[1];
rz(-3.1343757) q[2];
sx q[2];
rz(-1.3382148) q[2];
sx q[2];
rz(-0.54534334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73585549) q[1];
sx q[1];
rz(-0.89621004) q[1];
sx q[1];
rz(1.2184185) q[1];
rz(-pi) q[2];
rz(-2.7217322) q[3];
sx q[3];
rz(-1.6600939) q[3];
sx q[3];
rz(-1.4858703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-0.83690518) q[2];
sx q[2];
rz(-0.4591628) q[2];
rz(-3.0432213) q[3];
sx q[3];
rz(-1.8561074) q[3];
sx q[3];
rz(-1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(-0.14044811) q[0];
rz(-1.9239931) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(-1.5509031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2813148) q[0];
sx q[0];
rz(-1.4893805) q[0];
sx q[0];
rz(2.415262) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7920309) q[2];
sx q[2];
rz(-0.2436175) q[2];
sx q[2];
rz(-1.6175912) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.925736) q[1];
sx q[1];
rz(-1.3093776) q[1];
sx q[1];
rz(-2.5040313) q[1];
rz(-pi) q[2];
rz(0.012926558) q[3];
sx q[3];
rz(-2.2416246) q[3];
sx q[3];
rz(-2.304043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5156252) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(1.1930126) q[2];
rz(-0.48202816) q[3];
sx q[3];
rz(-0.41976443) q[3];
sx q[3];
rz(1.0198063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1358418) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(-0.76345481) q[1];
sx q[1];
rz(-1.8370942) q[1];
sx q[1];
rz(-2.7948517) q[1];
rz(-1.7557549) q[2];
sx q[2];
rz(-3.0228563) q[2];
sx q[2];
rz(-2.9517662) q[2];
rz(2.5478195) q[3];
sx q[3];
rz(-1.2111005) q[3];
sx q[3];
rz(2.0668798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
