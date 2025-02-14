OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5698009) q[0];
sx q[0];
rz(-2.2324012) q[0];
sx q[0];
rz(-1.5067014) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(6.7758898) q[1];
sx q[1];
rz(13.125782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5696379) q[0];
sx q[0];
rz(-0.35757726) q[0];
sx q[0];
rz(2.649369) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1790397) q[2];
sx q[2];
rz(-1.9843352) q[2];
sx q[2];
rz(-0.79038436) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8570429) q[1];
sx q[1];
rz(-0.27978292) q[1];
sx q[1];
rz(1.5955052) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2484965) q[3];
sx q[3];
rz(-0.71346766) q[3];
sx q[3];
rz(1.5763855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59014615) q[2];
sx q[2];
rz(-0.85942736) q[2];
sx q[2];
rz(1.0939416) q[2];
rz(1.9048196) q[3];
sx q[3];
rz(-1.9779132) q[3];
sx q[3];
rz(-2.8804603) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7673489) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(2.4236524) q[0];
rz(1.9473437) q[1];
sx q[1];
rz(-0.4078882) q[1];
sx q[1];
rz(-1.6494707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8555174) q[0];
sx q[0];
rz(-1.4465032) q[0];
sx q[0];
rz(0.13295023) q[0];
rz(-pi) q[1];
rz(-2.7020383) q[2];
sx q[2];
rz(-1.9157404) q[2];
sx q[2];
rz(2.3203691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7929505) q[1];
sx q[1];
rz(-1.1224282) q[1];
sx q[1];
rz(1.0690637) q[1];
rz(-pi) q[2];
rz(-0.82573311) q[3];
sx q[3];
rz(-1.7307708) q[3];
sx q[3];
rz(-0.79373327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0933928) q[2];
sx q[2];
rz(-2.3651249) q[2];
sx q[2];
rz(-0.36273599) q[2];
rz(-2.2705966) q[3];
sx q[3];
rz(-1.769915) q[3];
sx q[3];
rz(-1.3236275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094630346) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(0.56075019) q[0];
rz(2.1669855) q[1];
sx q[1];
rz(-0.86169306) q[1];
sx q[1];
rz(-2.9240756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1851674) q[0];
sx q[0];
rz(-2.5568058) q[0];
sx q[0];
rz(-0.034791481) q[0];
x q[1];
rz(3.0383695) q[2];
sx q[2];
rz(-2.0821314) q[2];
sx q[2];
rz(-2.5325506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88918658) q[1];
sx q[1];
rz(-1.4562618) q[1];
sx q[1];
rz(1.3088398) q[1];
x q[2];
rz(3.0514225) q[3];
sx q[3];
rz(-0.61093447) q[3];
sx q[3];
rz(-0.73874015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4983623) q[2];
sx q[2];
rz(-2.3537894) q[2];
sx q[2];
rz(-0.95334774) q[2];
rz(-1.0777473) q[3];
sx q[3];
rz(-1.6809623) q[3];
sx q[3];
rz(0.47422153) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.137602) q[0];
sx q[0];
rz(-2.9676262) q[0];
sx q[0];
rz(0.91651383) q[0];
rz(-0.42463955) q[1];
sx q[1];
rz(-1.3820796) q[1];
sx q[1];
rz(-2.0133846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4229168) q[0];
sx q[0];
rz(-0.87327086) q[0];
sx q[0];
rz(-0.6573702) q[0];
rz(-pi) q[1];
rz(-0.75707423) q[2];
sx q[2];
rz(-1.7221976) q[2];
sx q[2];
rz(0.30138563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27552106) q[1];
sx q[1];
rz(-2.0876608) q[1];
sx q[1];
rz(-1.9171153) q[1];
x q[2];
rz(-1.8966394) q[3];
sx q[3];
rz(-0.62281194) q[3];
sx q[3];
rz(2.4145221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26587036) q[2];
sx q[2];
rz(-1.6758726) q[2];
sx q[2];
rz(-1.9423368) q[2];
rz(1.1572329) q[3];
sx q[3];
rz(-2.3374989) q[3];
sx q[3];
rz(0.51586241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1542094) q[0];
sx q[0];
rz(-3.0976384) q[0];
sx q[0];
rz(-0.92535812) q[0];
rz(0.56272733) q[1];
sx q[1];
rz(-2.1268763) q[1];
sx q[1];
rz(-2.7664807) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5432388) q[0];
sx q[0];
rz(-1.028201) q[0];
sx q[0];
rz(-1.6981237) q[0];
x q[1];
rz(3.0882224) q[2];
sx q[2];
rz(-1.1359906) q[2];
sx q[2];
rz(-0.31440266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8598392) q[1];
sx q[1];
rz(-1.5956164) q[1];
sx q[1];
rz(2.9725705) q[1];
rz(-pi) q[2];
rz(-1.1697606) q[3];
sx q[3];
rz(-1.1479064) q[3];
sx q[3];
rz(-2.324375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.90317) q[2];
sx q[2];
rz(-0.91707245) q[2];
sx q[2];
rz(2.6440716) q[2];
rz(2.7072952) q[3];
sx q[3];
rz(-1.4562166) q[3];
sx q[3];
rz(-2.1527877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76181805) q[0];
sx q[0];
rz(-2.2803545) q[0];
sx q[0];
rz(0.27989835) q[0];
rz(-0.99413904) q[1];
sx q[1];
rz(-2.2810664) q[1];
sx q[1];
rz(2.5631189) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8461762) q[0];
sx q[0];
rz(-0.65447411) q[0];
sx q[0];
rz(-2.492948) q[0];
rz(2.4388588) q[2];
sx q[2];
rz(-1.8763855) q[2];
sx q[2];
rz(-1.0440799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.371767) q[1];
sx q[1];
rz(-2.3845362) q[1];
sx q[1];
rz(-0.28626059) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26831823) q[3];
sx q[3];
rz(-1.473241) q[3];
sx q[3];
rz(-2.8821111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0824739) q[2];
sx q[2];
rz(-2.5422577) q[2];
sx q[2];
rz(-2.2590051) q[2];
rz(2.5146218) q[3];
sx q[3];
rz(-2.4866703) q[3];
sx q[3];
rz(0.89491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249228) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(-0.00061568419) q[0];
rz(0.90283886) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(3.0965064) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49360291) q[0];
sx q[0];
rz(-1.6479074) q[0];
sx q[0];
rz(-0.48992975) q[0];
rz(-pi) q[1];
rz(0.74202027) q[2];
sx q[2];
rz(-1.4145523) q[2];
sx q[2];
rz(2.7458423) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0087456044) q[1];
sx q[1];
rz(-0.80423052) q[1];
sx q[1];
rz(-2.7102317) q[1];
x q[2];
rz(2.5302913) q[3];
sx q[3];
rz(-1.0795648) q[3];
sx q[3];
rz(-0.69913759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8467329) q[2];
sx q[2];
rz(-0.555987) q[2];
sx q[2];
rz(-0.9168469) q[2];
rz(2.2461241) q[3];
sx q[3];
rz(-1.5119036) q[3];
sx q[3];
rz(1.2112613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13670066) q[0];
sx q[0];
rz(-0.08575511) q[0];
sx q[0];
rz(-2.6766747) q[0];
rz(1.1888986) q[1];
sx q[1];
rz(-1.7949972) q[1];
sx q[1];
rz(2.7704923) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2087323) q[0];
sx q[0];
rz(-1.1602872) q[0];
sx q[0];
rz(-2.6952144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6109054) q[2];
sx q[2];
rz(-0.75119392) q[2];
sx q[2];
rz(-1.4152648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6976561) q[1];
sx q[1];
rz(-2.3015226) q[1];
sx q[1];
rz(-0.96644641) q[1];
rz(0.38406541) q[3];
sx q[3];
rz(-1.8397985) q[3];
sx q[3];
rz(-0.88685461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4798639) q[2];
sx q[2];
rz(-2.3532545) q[2];
sx q[2];
rz(-2.7313477) q[2];
rz(1.6346301) q[3];
sx q[3];
rz(-2.7362636) q[3];
sx q[3];
rz(-3.098367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048454408) q[0];
sx q[0];
rz(-0.52670902) q[0];
sx q[0];
rz(-0.17573892) q[0];
rz(-2.3161092) q[1];
sx q[1];
rz(-1.261542) q[1];
sx q[1];
rz(2.3186191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776445) q[0];
sx q[0];
rz(-0.86714449) q[0];
sx q[0];
rz(-2.4022341) q[0];
rz(-1.9972598) q[2];
sx q[2];
rz(-1.8606295) q[2];
sx q[2];
rz(-2.5215931) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8205679) q[1];
sx q[1];
rz(-1.9181644) q[1];
sx q[1];
rz(-2.419761) q[1];
rz(-pi) q[2];
rz(-1.6552583) q[3];
sx q[3];
rz(-0.41197398) q[3];
sx q[3];
rz(0.99124817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2814111) q[2];
sx q[2];
rz(-1.1450291) q[2];
sx q[2];
rz(-0.54110503) q[2];
rz(-2.7152854) q[3];
sx q[3];
rz(-0.46190327) q[3];
sx q[3];
rz(0.84686744) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68778872) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(2.8059106) q[0];
rz(3.1188534) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(2.4679599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84647932) q[0];
sx q[0];
rz(-1.9178977) q[0];
sx q[0];
rz(-0.62008574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9433374) q[2];
sx q[2];
rz(-0.47706826) q[2];
sx q[2];
rz(-1.4587634) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9424098) q[1];
sx q[1];
rz(-1.3681776) q[1];
sx q[1];
rz(2.1658648) q[1];
rz(-pi) q[2];
rz(1.4309747) q[3];
sx q[3];
rz(-1.9019902) q[3];
sx q[3];
rz(-3.0013623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.444904) q[2];
sx q[2];
rz(-2.1477369) q[2];
sx q[2];
rz(-2.135684) q[2];
rz(-0.78399793) q[3];
sx q[3];
rz(-2.8204212) q[3];
sx q[3];
rz(-1.706749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0189331) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(-2.7574273) q[1];
sx q[1];
rz(-1.1687678) q[1];
sx q[1];
rz(1.6484177) q[1];
rz(-0.76569414) q[2];
sx q[2];
rz(-0.29065825) q[2];
sx q[2];
rz(-2.7343168) q[2];
rz(2.7586423) q[3];
sx q[3];
rz(-0.89385116) q[3];
sx q[3];
rz(1.6123621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
