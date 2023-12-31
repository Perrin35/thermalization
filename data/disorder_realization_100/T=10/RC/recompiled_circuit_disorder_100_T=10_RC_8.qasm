OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(1.4555414) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769387) q[0];
sx q[0];
rz(-1.1507478) q[0];
sx q[0];
rz(2.7647892) q[0];
rz(-0.86906616) q[2];
sx q[2];
rz(-1.367374) q[2];
sx q[2];
rz(-2.2609401) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9283596) q[1];
sx q[1];
rz(-1.8395437) q[1];
sx q[1];
rz(2.1409722) q[1];
rz(-pi) q[2];
rz(-3.0936712) q[3];
sx q[3];
rz(-1.5600292) q[3];
sx q[3];
rz(-2.3445232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(-0.39164266) q[2];
rz(-3.0018905) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.93719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6635839) q[0];
sx q[0];
rz(-1.9933797) q[0];
sx q[0];
rz(2.708486) q[0];
x q[1];
rz(1.3538829) q[2];
sx q[2];
rz(-0.7173983) q[2];
sx q[2];
rz(2.4153828) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.345563) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(1.4644535) q[1];
rz(0.28537206) q[3];
sx q[3];
rz(-0.91802363) q[3];
sx q[3];
rz(-2.5170346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.5725117) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(-0.61808008) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.9504257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55710775) q[0];
sx q[0];
rz(-2.7911721) q[0];
sx q[0];
rz(-1.114203) q[0];
rz(-2.4345258) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(2.6235839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9900073) q[1];
sx q[1];
rz(-1.7813492) q[1];
sx q[1];
rz(-1.3596119) q[1];
rz(2.5997945) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1742192) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8106666) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(1.7640242) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-0.92299443) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3337293) q[0];
sx q[0];
rz(-2.6834052) q[0];
sx q[0];
rz(-2.3966167) q[0];
rz(-pi) q[1];
x q[1];
rz(0.043945233) q[2];
sx q[2];
rz(-1.5825795) q[2];
sx q[2];
rz(0.95566434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.1715602) q[1];
rz(-0.8664341) q[3];
sx q[3];
rz(-1.5537062) q[3];
sx q[3];
rz(-0.41364241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4541645) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-0.50951177) q[2];
rz(-2.4984958) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(2.2221785) q[0];
rz(-2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.3607508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0589941) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(-pi) q[1];
rz(2.1941357) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(-2.5156977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.67140019) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(2.7027674) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0642062) q[3];
sx q[3];
rz(-0.97273982) q[3];
sx q[3];
rz(3.0706881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.8615287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(0.70760977) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(2.8663666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.595364) q[0];
sx q[0];
rz(-1.9018136) q[0];
sx q[0];
rz(2.0966895) q[0];
rz(2.2116824) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(-0.90604679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3616398) q[1];
sx q[1];
rz(-1.0859023) q[1];
sx q[1];
rz(1.9655025) q[1];
x q[2];
rz(2.1915057) q[3];
sx q[3];
rz(-1.5035149) q[3];
sx q[3];
rz(2.7954807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-0.71969676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450711) q[0];
sx q[0];
rz(-0.70598999) q[0];
sx q[0];
rz(2.6939874) q[0];
rz(-pi) q[1];
rz(0.72549707) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(3.0692435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.91832671) q[1];
sx q[1];
rz(-2.1152788) q[1];
sx q[1];
rz(1.6225999) q[1];
x q[2];
rz(1.3004488) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(-2.8519252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197501) q[0];
sx q[0];
rz(-1.9103721) q[0];
sx q[0];
rz(1.3440488) q[0];
rz(2.9096793) q[2];
sx q[2];
rz(-0.70966087) q[2];
sx q[2];
rz(-2.2556925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9547792) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(0.38444744) q[1];
rz(-0.90663119) q[3];
sx q[3];
rz(-2.2615221) q[3];
sx q[3];
rz(3.0674792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(1.2773369) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(0.87798464) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.3605798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.774705) q[0];
sx q[0];
rz(-2.5686712) q[0];
sx q[0];
rz(1.5013055) q[0];
rz(-0.94366818) q[2];
sx q[2];
rz(-1.6086279) q[2];
sx q[2];
rz(-0.19247069) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18371836) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(-2.3493489) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9079886) q[3];
sx q[3];
rz(-2.7003151) q[3];
sx q[3];
rz(0.518276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.517841) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(3.1023074) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601446) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(-3.0597276) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1245329) q[2];
sx q[2];
rz(-2.8327201) q[2];
sx q[2];
rz(2.942254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7457909) q[1];
sx q[1];
rz(-1.5862641) q[1];
sx q[1];
rz(-1.323911) q[1];
x q[2];
rz(-2.2640957) q[3];
sx q[3];
rz(-0.96713669) q[3];
sx q[3];
rz(0.391215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7211192) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532886) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(0.2164671) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
rz(1.5218232) q[3];
sx q[3];
rz(-2.5544142) q[3];
sx q[3];
rz(2.1094473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
