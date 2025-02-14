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
rz(0.94035971) q[0];
sx q[0];
rz(-2.0191963) q[0];
sx q[0];
rz(2.9474592) q[0];
rz(1.0951207) q[1];
sx q[1];
rz(-1.6835338) q[1];
sx q[1];
rz(2.3966052) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4177822) q[0];
sx q[0];
rz(-3.0273154) q[0];
sx q[0];
rz(2.6835915) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3294747) q[2];
sx q[2];
rz(-1.0924763) q[2];
sx q[2];
rz(1.1253192) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2304334) q[1];
sx q[1];
rz(-1.8130366) q[1];
sx q[1];
rz(-1.6723067) q[1];
rz(0.93985112) q[3];
sx q[3];
rz(-2.5999024) q[3];
sx q[3];
rz(1.9909137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8038883) q[2];
sx q[2];
rz(-1.3362198) q[2];
sx q[2];
rz(1.3275006) q[2];
rz(1.7727857) q[3];
sx q[3];
rz(-2.8076706) q[3];
sx q[3];
rz(-0.22148618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1526445) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(-0.055135559) q[0];
rz(-0.70471835) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(-1.0447186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7009788) q[0];
sx q[0];
rz(-1.199628) q[0];
sx q[0];
rz(1.2775902) q[0];
x q[1];
rz(-0.98590322) q[2];
sx q[2];
rz(-1.0828746) q[2];
sx q[2];
rz(0.019854644) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4028496) q[1];
sx q[1];
rz(-0.84058769) q[1];
sx q[1];
rz(0.56096803) q[1];
x q[2];
rz(3.0468076) q[3];
sx q[3];
rz(-1.6322502) q[3];
sx q[3];
rz(0.8549785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3220871) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(0.78754887) q[3];
sx q[3];
rz(-1.4044263) q[3];
sx q[3];
rz(-2.0866086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2128485) q[0];
sx q[0];
rz(-2.4725914) q[0];
sx q[0];
rz(-2.4231732) q[0];
rz(0.33572117) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(-2.1885923) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7756726) q[0];
sx q[0];
rz(-1.6073415) q[0];
sx q[0];
rz(2.1120872) q[0];
x q[1];
rz(1.3173872) q[2];
sx q[2];
rz(-2.6077787) q[2];
sx q[2];
rz(1.6718503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53297199) q[1];
sx q[1];
rz(-2.2257881) q[1];
sx q[1];
rz(2.6650732) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5312214) q[3];
sx q[3];
rz(-2.000375) q[3];
sx q[3];
rz(-2.3342032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(2.4887264) q[2];
rz(-0.79814923) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(-0.64364141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2586486) q[0];
sx q[0];
rz(-0.84027165) q[0];
sx q[0];
rz(0.81664455) q[0];
rz(0.69149292) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(-0.74849558) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25835055) q[0];
sx q[0];
rz(-1.878698) q[0];
sx q[0];
rz(1.3624205) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9709936) q[2];
sx q[2];
rz(-0.23175254) q[2];
sx q[2];
rz(-2.0830926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4481727) q[1];
sx q[1];
rz(-2.2058371) q[1];
sx q[1];
rz(0.5902115) q[1];
x q[2];
rz(0.17552142) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(2.1057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3889435) q[2];
sx q[2];
rz(-2.3708673) q[2];
sx q[2];
rz(-2.3469817) q[2];
rz(1.7157582) q[3];
sx q[3];
rz(-1.040193) q[3];
sx q[3];
rz(-0.37253255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(2.6439731) q[0];
rz(-0.76958641) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(2.41113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32942922) q[0];
sx q[0];
rz(-2.4162628) q[0];
sx q[0];
rz(2.8149288) q[0];
rz(-pi) q[1];
rz(2.8874997) q[2];
sx q[2];
rz(-0.95839082) q[2];
sx q[2];
rz(-1.8191561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9047315) q[1];
sx q[1];
rz(-1.3610468) q[1];
sx q[1];
rz(2.3244203) q[1];
rz(-pi) q[2];
rz(-0.83010556) q[3];
sx q[3];
rz(-2.5200994) q[3];
sx q[3];
rz(1.1661539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89255014) q[2];
sx q[2];
rz(-1.6232792) q[2];
sx q[2];
rz(0.45336938) q[2];
rz(1.8306277) q[3];
sx q[3];
rz(-0.1736719) q[3];
sx q[3];
rz(-3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740358) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(2.0454316) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(-0.62166628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8209648) q[0];
sx q[0];
rz(-0.85658565) q[0];
sx q[0];
rz(-1.5498534) q[0];
x q[1];
rz(-1.4716343) q[2];
sx q[2];
rz(-1.9648203) q[2];
sx q[2];
rz(-2.9539915) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0677133) q[1];
sx q[1];
rz(-0.82649684) q[1];
sx q[1];
rz(-1.9113051) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4765763) q[3];
sx q[3];
rz(-1.1248853) q[3];
sx q[3];
rz(-2.8578106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93214503) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(-1.4383593) q[2];
rz(2.1608593) q[3];
sx q[3];
rz(-1.2659729) q[3];
sx q[3];
rz(1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.257523) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(-0.32514611) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-2.1279361) q[1];
sx q[1];
rz(1.315717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60916736) q[0];
sx q[0];
rz(-2.2497674) q[0];
sx q[0];
rz(0.15773134) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1431214) q[2];
sx q[2];
rz(-1.1954167) q[2];
sx q[2];
rz(-1.4294415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.084242227) q[1];
sx q[1];
rz(-0.92613832) q[1];
sx q[1];
rz(-2.0688384) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93856424) q[3];
sx q[3];
rz(-1.5364416) q[3];
sx q[3];
rz(-1.1704579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.331984) q[2];
sx q[2];
rz(-2.1029682) q[2];
sx q[2];
rz(1.5254376) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.6922035) q[3];
sx q[3];
rz(-1.6239032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759393) q[0];
sx q[0];
rz(-2.0015367) q[0];
sx q[0];
rz(-0.64919382) q[0];
rz(2.3340885) q[1];
sx q[1];
rz(-1.1367926) q[1];
sx q[1];
rz(-1.753122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8313468) q[0];
sx q[0];
rz(-2.5355336) q[0];
sx q[0];
rz(1.640727) q[0];
rz(-pi) q[1];
rz(-1.1325324) q[2];
sx q[2];
rz(-0.31719855) q[2];
sx q[2];
rz(0.68480051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5538881) q[1];
sx q[1];
rz(-0.18408751) q[1];
sx q[1];
rz(1.9485997) q[1];
x q[2];
rz(-1.3598594) q[3];
sx q[3];
rz(-1.9550793) q[3];
sx q[3];
rz(-2.3944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2041152) q[2];
sx q[2];
rz(-2.9948339) q[2];
sx q[2];
rz(0.80782962) q[2];
rz(-2.4701123) q[3];
sx q[3];
rz(-1.6356133) q[3];
sx q[3];
rz(-1.5427264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.2272187) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(-2.6887023) q[0];
rz(0.29531404) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(-1.7426851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5495396) q[0];
sx q[0];
rz(-2.4045134) q[0];
sx q[0];
rz(-1.6960742) q[0];
x q[1];
rz(-2.9851297) q[2];
sx q[2];
rz(-1.0623806) q[2];
sx q[2];
rz(3.1389295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4537482) q[1];
sx q[1];
rz(-1.1198938) q[1];
sx q[1];
rz(-1.1503673) q[1];
rz(-pi) q[2];
rz(-1.6803375) q[3];
sx q[3];
rz(-2.9432812) q[3];
sx q[3];
rz(1.2480145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7704775) q[2];
sx q[2];
rz(-2.2687843) q[2];
sx q[2];
rz(2.5650043) q[2];
rz(2.9234486) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(-3.0531378) q[0];
rz(0.5312008) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(-1.521135) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59014295) q[0];
sx q[0];
rz(-1.1781173) q[0];
sx q[0];
rz(3.0513502) q[0];
rz(-pi) q[1];
rz(3.0938153) q[2];
sx q[2];
rz(-2.0476404) q[2];
sx q[2];
rz(0.58281118) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.759821) q[1];
sx q[1];
rz(-2.9050937) q[1];
sx q[1];
rz(2.8245154) q[1];
rz(-pi) q[2];
rz(1.985667) q[3];
sx q[3];
rz(-2.1262827) q[3];
sx q[3];
rz(-1.5145356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8497808) q[2];
sx q[2];
rz(-1.7917289) q[2];
sx q[2];
rz(-1.6144217) q[2];
rz(-1.269086) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.108718) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(-2.6678008) q[1];
sx q[1];
rz(-1.5693249) q[1];
sx q[1];
rz(1.5798461) q[1];
rz(-2.5529591) q[2];
sx q[2];
rz(-0.6856754) q[2];
sx q[2];
rz(-2.6720809) q[2];
rz(-0.56415767) q[3];
sx q[3];
rz(-2.0750792) q[3];
sx q[3];
rz(0.32479494) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
