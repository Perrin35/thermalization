OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44242087) q[0];
sx q[0];
rz(-2.3306263) q[0];
sx q[0];
rz(2.6851658) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(2.9589597) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62902495) q[0];
sx q[0];
rz(-2.1898666) q[0];
sx q[0];
rz(-2.7946212) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86948189) q[2];
sx q[2];
rz(-1.9449781) q[2];
sx q[2];
rz(2.2187198) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7108954) q[1];
sx q[1];
rz(-0.17696807) q[1];
sx q[1];
rz(-2.9269993) q[1];
x q[2];
rz(-1.9630646) q[3];
sx q[3];
rz(-1.0844106) q[3];
sx q[3];
rz(1.148759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7543588) q[2];
sx q[2];
rz(-2.0626455) q[2];
sx q[2];
rz(0.31164247) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-2.8897372) q[3];
sx q[3];
rz(-0.18865147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(1.9022994) q[0];
rz(2.7481825) q[1];
sx q[1];
rz(-2.0937803) q[1];
sx q[1];
rz(-2.1107103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9728972) q[0];
sx q[0];
rz(-2.4542232) q[0];
sx q[0];
rz(1.0052135) q[0];
rz(-0.26427697) q[2];
sx q[2];
rz(-1.7925996) q[2];
sx q[2];
rz(1.4724195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6759265) q[1];
sx q[1];
rz(-0.74518004) q[1];
sx q[1];
rz(1.0449157) q[1];
rz(-2.2575284) q[3];
sx q[3];
rz(-0.86694781) q[3];
sx q[3];
rz(0.91020179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82765141) q[2];
sx q[2];
rz(-0.88259077) q[2];
sx q[2];
rz(-2.7871056) q[2];
rz(-0.83141023) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(-1.6262511) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77446929) q[0];
sx q[0];
rz(-0.18018436) q[0];
sx q[0];
rz(-0.48429504) q[0];
rz(0.88227415) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(-0.54214111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4216293) q[0];
sx q[0];
rz(-2.6609592) q[0];
sx q[0];
rz(-2.292657) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32140478) q[2];
sx q[2];
rz(-1.5084477) q[2];
sx q[2];
rz(-2.7624109) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.650108) q[1];
sx q[1];
rz(-1.3940147) q[1];
sx q[1];
rz(-0.80727908) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8081424) q[3];
sx q[3];
rz(-1.508731) q[3];
sx q[3];
rz(-0.82236034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.813628) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(-2.1072809) q[2];
rz(1.7799001) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(2.5907607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4807602) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(-2.0939636) q[0];
rz(1.818559) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(-0.38527647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40257257) q[0];
sx q[0];
rz(-2.077335) q[0];
sx q[0];
rz(-2.6468524) q[0];
rz(-0.71661894) q[2];
sx q[2];
rz(-1.5138549) q[2];
sx q[2];
rz(-0.54887923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1487354) q[1];
sx q[1];
rz(-2.9210733) q[1];
sx q[1];
rz(0.33111568) q[1];
rz(-pi) q[2];
rz(3.0317467) q[3];
sx q[3];
rz(-0.89021909) q[3];
sx q[3];
rz(0.79119134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72383991) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(0.27457944) q[2];
rz(2.5802021) q[3];
sx q[3];
rz(-1.0654819) q[3];
sx q[3];
rz(-1.6339462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0059589) q[0];
sx q[0];
rz(-2.5649286) q[0];
sx q[0];
rz(0.86719257) q[0];
rz(2.7658956) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(1.2295178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4081683) q[0];
sx q[0];
rz(-0.46981341) q[0];
sx q[0];
rz(-2.7759659) q[0];
x q[1];
rz(2.0784573) q[2];
sx q[2];
rz(-2.5637321) q[2];
sx q[2];
rz(0.40921989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2676437) q[1];
sx q[1];
rz(-0.76890828) q[1];
sx q[1];
rz(2.257009) q[1];
rz(-pi) q[2];
rz(2.8035012) q[3];
sx q[3];
rz(-1.1792569) q[3];
sx q[3];
rz(0.028926802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.025297252) q[2];
sx q[2];
rz(-1.2325492) q[2];
sx q[2];
rz(0.06289014) q[2];
rz(-0.40999117) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(-2.9819152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9922239) q[0];
sx q[0];
rz(-3.0719482) q[0];
sx q[0];
rz(-0.83576354) q[0];
rz(1.958485) q[1];
sx q[1];
rz(-1.2703905) q[1];
sx q[1];
rz(0.74367181) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8643657) q[0];
sx q[0];
rz(-1.0334618) q[0];
sx q[0];
rz(-1.9242084) q[0];
rz(-1.7149107) q[2];
sx q[2];
rz(-1.5169355) q[2];
sx q[2];
rz(0.28648057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89827713) q[1];
sx q[1];
rz(-2.26663) q[1];
sx q[1];
rz(-0.23855539) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25571574) q[3];
sx q[3];
rz(-1.6220495) q[3];
sx q[3];
rz(1.071614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30367294) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(0.79353235) q[2];
rz(-2.7225336) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(-0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917052) q[0];
sx q[0];
rz(-2.1623623) q[0];
sx q[0];
rz(-1.9784084) q[0];
rz(-2.361182) q[1];
sx q[1];
rz(-2.8051832) q[1];
sx q[1];
rz(-1.6132678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0931143) q[0];
sx q[0];
rz(-2.265904) q[0];
sx q[0];
rz(3.0570875) q[0];
rz(-pi) q[1];
rz(2.9588671) q[2];
sx q[2];
rz(-1.2264681) q[2];
sx q[2];
rz(-1.2862213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1133729) q[1];
sx q[1];
rz(-1.372678) q[1];
sx q[1];
rz(0.021620167) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9933543) q[3];
sx q[3];
rz(-1.3179639) q[3];
sx q[3];
rz(-1.8971083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0344326) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(3.072928) q[2];
rz(2.5717403) q[3];
sx q[3];
rz(-2.6611501) q[3];
sx q[3];
rz(-0.3705875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625075) q[0];
sx q[0];
rz(-0.67254368) q[0];
sx q[0];
rz(-1.3154718) q[0];
rz(1.8585809) q[1];
sx q[1];
rz(-0.43671572) q[1];
sx q[1];
rz(-3.0924996) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6460719) q[0];
sx q[0];
rz(-1.4804513) q[0];
sx q[0];
rz(3.0642088) q[0];
rz(-pi) q[1];
rz(0.53907271) q[2];
sx q[2];
rz(-2.1207715) q[2];
sx q[2];
rz(1.0895188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1421353) q[1];
sx q[1];
rz(-1.9519567) q[1];
sx q[1];
rz(0.8183523) q[1];
rz(3.0144948) q[3];
sx q[3];
rz(-0.93080257) q[3];
sx q[3];
rz(0.52332969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38356885) q[2];
sx q[2];
rz(-0.87791666) q[2];
sx q[2];
rz(0.54086584) q[2];
rz(2.0567549) q[3];
sx q[3];
rz(-0.70298755) q[3];
sx q[3];
rz(-1.7613523) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.908602) q[0];
sx q[0];
rz(-0.99642307) q[0];
sx q[0];
rz(-0.10502271) q[0];
rz(0.54166334) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(-2.7856316) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0070156688) q[0];
sx q[0];
rz(-1.4316214) q[0];
sx q[0];
rz(-0.037875847) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84635205) q[2];
sx q[2];
rz(-1.0484107) q[2];
sx q[2];
rz(1.9047996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7193422) q[1];
sx q[1];
rz(-1.0021035) q[1];
sx q[1];
rz(2.689792) q[1];
rz(-1.2508873) q[3];
sx q[3];
rz(-1.6875132) q[3];
sx q[3];
rz(-3.0206783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2124704) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(-1.8222202) q[2];
rz(1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(-2.1738079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(1.6784278) q[0];
rz(-1.6819008) q[1];
sx q[1];
rz(-1.1594783) q[1];
sx q[1];
rz(2.7957338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1347156) q[0];
sx q[0];
rz(-1.6505989) q[0];
sx q[0];
rz(0.047264506) q[0];
rz(-1.6357291) q[2];
sx q[2];
rz(-1.254515) q[2];
sx q[2];
rz(2.8313178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74259669) q[1];
sx q[1];
rz(-2.1861939) q[1];
sx q[1];
rz(1.0091429) q[1];
rz(-0.34158753) q[3];
sx q[3];
rz(-0.95529592) q[3];
sx q[3];
rz(-3.0581491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.8701376) q[2];
sx q[2];
rz(2.8806809) q[2];
rz(1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3093001) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(-0.55083864) q[1];
sx q[1];
rz(-1.486634) q[1];
sx q[1];
rz(-2.6015729) q[1];
rz(1.3030686) q[2];
sx q[2];
rz(-2.3961551) q[2];
sx q[2];
rz(1.2533631) q[2];
rz(2.8367219) q[3];
sx q[3];
rz(-2.3875227) q[3];
sx q[3];
rz(1.0878121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
