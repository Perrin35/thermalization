OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.85668808) q[0];
sx q[0];
rz(-1.3311128) q[0];
sx q[0];
rz(-2.3321505) q[0];
rz(0.59960214) q[1];
sx q[1];
rz(-2.0166346) q[1];
sx q[1];
rz(-0.52176276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7441352) q[0];
sx q[0];
rz(-1.43072) q[0];
sx q[0];
rz(0.0011573275) q[0];
rz(-pi) q[1];
rz(-2.6894301) q[2];
sx q[2];
rz(-0.60930646) q[2];
sx q[2];
rz(0.094560187) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45568902) q[1];
sx q[1];
rz(-1.1215116) q[1];
sx q[1];
rz(-0.30251578) q[1];
rz(-2.605569) q[3];
sx q[3];
rz(-1.496435) q[3];
sx q[3];
rz(-2.1002567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1521505) q[2];
sx q[2];
rz(-2.4336954) q[2];
sx q[2];
rz(1.36261) q[2];
rz(1.4710434) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(-2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489557) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(1.0697399) q[0];
rz(2.3906129) q[1];
sx q[1];
rz(-0.90194482) q[1];
sx q[1];
rz(0.78838563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716924) q[0];
sx q[0];
rz(-2.4909752) q[0];
sx q[0];
rz(0.31221892) q[0];
x q[1];
rz(-0.66547243) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(1.9201345) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3020116) q[1];
sx q[1];
rz(-2.5835648) q[1];
sx q[1];
rz(1.9424979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51075824) q[3];
sx q[3];
rz(-1.184297) q[3];
sx q[3];
rz(-3.136477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.02701935) q[2];
sx q[2];
rz(-2.7832649) q[2];
sx q[2];
rz(-0.9066073) q[2];
rz(-1.00057) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(-0.51893836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(-1.3038127) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(1.7283641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50158635) q[0];
sx q[0];
rz(-1.4297856) q[0];
sx q[0];
rz(1.4835351) q[0];
rz(-pi) q[1];
rz(-1.1151074) q[2];
sx q[2];
rz(-2.1564061) q[2];
sx q[2];
rz(-0.93916303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6709137) q[1];
sx q[1];
rz(-1.7470198) q[1];
sx q[1];
rz(0.46609585) q[1];
x q[2];
rz(1.71536) q[3];
sx q[3];
rz(-0.72025296) q[3];
sx q[3];
rz(1.4831869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.383519) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-3.0124532) q[2];
rz(0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(-2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(-2.2953798) q[0];
rz(0.57202488) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(-1.2342854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0219824) q[0];
sx q[0];
rz(-1.9024769) q[0];
sx q[0];
rz(-1.4681547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1051142) q[2];
sx q[2];
rz(-0.52157516) q[2];
sx q[2];
rz(1.2273481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14951359) q[1];
sx q[1];
rz(-1.7798335) q[1];
sx q[1];
rz(1.4884218) q[1];
x q[2];
rz(2.3448496) q[3];
sx q[3];
rz(-1.1827743) q[3];
sx q[3];
rz(1.6643559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7531551) q[2];
sx q[2];
rz(-0.93418241) q[2];
sx q[2];
rz(1.4446806) q[2];
rz(1.2706903) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020849) q[0];
sx q[0];
rz(-1.6725699) q[0];
sx q[0];
rz(2.0462659) q[0];
rz(0.0630088) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(0.94620401) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444855) q[0];
sx q[0];
rz(-1.5427663) q[0];
sx q[0];
rz(2.1155684) q[0];
rz(-pi) q[1];
rz(2.6465552) q[2];
sx q[2];
rz(-2.4433141) q[2];
sx q[2];
rz(-0.78294047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9621219) q[1];
sx q[1];
rz(-1.0596766) q[1];
sx q[1];
rz(1.5764135) q[1];
rz(1.0466411) q[3];
sx q[3];
rz(-0.44455179) q[3];
sx q[3];
rz(-1.1994282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63518628) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(1.6335454) q[2];
rz(0.15277282) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(-0.93301409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372811) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(-3.0142504) q[1];
sx q[1];
rz(-1.2510977) q[1];
sx q[1];
rz(0.88643518) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3484257) q[0];
sx q[0];
rz(-2.0901449) q[0];
sx q[0];
rz(-1.5898806) q[0];
x q[1];
rz(1.7753052) q[2];
sx q[2];
rz(-1.5123774) q[2];
sx q[2];
rz(1.8522592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97798367) q[1];
sx q[1];
rz(-2.915288) q[1];
sx q[1];
rz(-1.3127021) q[1];
x q[2];
rz(1.4599019) q[3];
sx q[3];
rz(-1.4023047) q[3];
sx q[3];
rz(3.0805825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2966557) q[2];
sx q[2];
rz(-2.3462494) q[2];
sx q[2];
rz(-0.33357683) q[2];
rz(0.9345471) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(-2.2540588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6880671) q[0];
sx q[0];
rz(-1.4997361) q[0];
sx q[0];
rz(-2.8330579) q[0];
rz(-0.60641369) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(-2.6920614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21643695) q[0];
sx q[0];
rz(-2.0846268) q[0];
sx q[0];
rz(-1.6440744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1029195) q[2];
sx q[2];
rz(-0.84005594) q[2];
sx q[2];
rz(2.8401224) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7923741) q[1];
sx q[1];
rz(-1.6728396) q[1];
sx q[1];
rz(-0.95178638) q[1];
rz(0.94664871) q[3];
sx q[3];
rz(-0.80393857) q[3];
sx q[3];
rz(-2.338666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0799847) q[2];
sx q[2];
rz(-1.6500762) q[2];
sx q[2];
rz(-1.9645346) q[2];
rz(-0.70147771) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(-2.285932) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45706448) q[0];
sx q[0];
rz(-0.51790154) q[0];
sx q[0];
rz(1.49217) q[0];
rz(2.0860784) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(-0.42748705) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9857835) q[0];
sx q[0];
rz(-2.0844578) q[0];
sx q[0];
rz(2.2751121) q[0];
rz(-pi) q[1];
rz(-1.7321109) q[2];
sx q[2];
rz(-1.1448793) q[2];
sx q[2];
rz(-2.6420455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0574574) q[1];
sx q[1];
rz(-1.5930719) q[1];
sx q[1];
rz(-0.31288625) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70722039) q[3];
sx q[3];
rz(-0.40509352) q[3];
sx q[3];
rz(1.1680548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18275729) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(-1.9367564) q[2];
rz(-0.13293535) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8590915) q[0];
sx q[0];
rz(-1.7318672) q[0];
sx q[0];
rz(-2.6891563) q[0];
rz(-0.85894194) q[1];
sx q[1];
rz(-2.5622538) q[1];
sx q[1];
rz(-1.4525684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4005139) q[0];
sx q[0];
rz(-0.79325926) q[0];
sx q[0];
rz(0.69852065) q[0];
rz(-0.33272393) q[2];
sx q[2];
rz(-2.1004268) q[2];
sx q[2];
rz(-0.67654787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0060785) q[1];
sx q[1];
rz(-1.0648512) q[1];
sx q[1];
rz(-3.1385026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9296367) q[3];
sx q[3];
rz(-1.9178041) q[3];
sx q[3];
rz(1.8458091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9026044) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(-0.39546173) q[2];
rz(-2.0579193) q[3];
sx q[3];
rz(-2.5189221) q[3];
sx q[3];
rz(-1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.285242) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(1.850199) q[0];
rz(-0.77766386) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(-1.2819598) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4829452) q[0];
sx q[0];
rz(-1.3532257) q[0];
sx q[0];
rz(-1.3402142) q[0];
rz(-pi) q[1];
rz(0.28996946) q[2];
sx q[2];
rz(-0.66807884) q[2];
sx q[2];
rz(-1.5529322) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3471372) q[1];
sx q[1];
rz(-1.0442808) q[1];
sx q[1];
rz(-1.8361676) q[1];
x q[2];
rz(-2.090023) q[3];
sx q[3];
rz(-1.4759721) q[3];
sx q[3];
rz(-0.42025305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97734863) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(-0.39883167) q[2];
rz(0.93419689) q[3];
sx q[3];
rz(-0.7312921) q[3];
sx q[3];
rz(-0.092223316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4451404) q[0];
sx q[0];
rz(-0.91309375) q[0];
sx q[0];
rz(0.0050807411) q[0];
rz(-2.483881) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(-1.896454) q[2];
sx q[2];
rz(-2.3477868) q[2];
sx q[2];
rz(-3.1127047) q[2];
rz(-2.1565831) q[3];
sx q[3];
rz(-1.6228709) q[3];
sx q[3];
rz(-3.090938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
