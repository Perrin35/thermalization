OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23842731) q[0];
sx q[0];
rz(-2.981346) q[0];
sx q[0];
rz(1.2877553) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(2.3828659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903126) q[0];
sx q[0];
rz(-2.7669124) q[0];
sx q[0];
rz(2.4640342) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1333732) q[2];
sx q[2];
rz(-2.6133399) q[2];
sx q[2];
rz(-0.6483486) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.030144) q[1];
sx q[1];
rz(-1.2093102) q[1];
sx q[1];
rz(0.8059146) q[1];
x q[2];
rz(-0.82697273) q[3];
sx q[3];
rz(-2.0872413) q[3];
sx q[3];
rz(1.4475105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(2.63499) q[2];
rz(1.5233585) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(-3.0860331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697206) q[0];
sx q[0];
rz(-0.46590081) q[0];
sx q[0];
rz(0.063657612) q[0];
rz(1.1977389) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(-1.0187842) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5757345) q[0];
sx q[0];
rz(-1.3745148) q[0];
sx q[0];
rz(1.7193166) q[0];
rz(1.3345785) q[2];
sx q[2];
rz(-1.2389606) q[2];
sx q[2];
rz(-0.97768367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19551046) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(2.4184188) q[1];
rz(3.0080454) q[3];
sx q[3];
rz(-0.822328) q[3];
sx q[3];
rz(-0.68822569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7055052) q[2];
sx q[2];
rz(-2.5477396) q[2];
sx q[2];
rz(0.84211055) q[2];
rz(1.0574794) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975526) q[0];
sx q[0];
rz(-1.0201447) q[0];
sx q[0];
rz(1.9496339) q[0];
rz(0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.8398197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98604964) q[0];
sx q[0];
rz(-1.0240004) q[0];
sx q[0];
rz(0.56122551) q[0];
x q[1];
rz(-2.9495732) q[2];
sx q[2];
rz(-2.9387967) q[2];
sx q[2];
rz(-1.6209728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5980581) q[1];
sx q[1];
rz(-0.50331668) q[1];
sx q[1];
rz(-0.37737585) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93940063) q[3];
sx q[3];
rz(-2.7223848) q[3];
sx q[3];
rz(1.7260176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7630345) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(-3.0540826) q[2];
rz(-1.3148426) q[3];
sx q[3];
rz(-2.0851236) q[3];
sx q[3];
rz(0.99353138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4104376) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(-0.80899578) q[0];
rz(1.0357098) q[1];
sx q[1];
rz(-0.46329841) q[1];
sx q[1];
rz(2.1083924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3147723) q[0];
sx q[0];
rz(-2.2345532) q[0];
sx q[0];
rz(-0.04962595) q[0];
x q[1];
rz(1.0339917) q[2];
sx q[2];
rz(-1.9516757) q[2];
sx q[2];
rz(1.8465259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3174008) q[1];
sx q[1];
rz(-0.5017952) q[1];
sx q[1];
rz(0.3806033) q[1];
rz(-0.54368005) q[3];
sx q[3];
rz(-1.1515695) q[3];
sx q[3];
rz(-0.70862428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0741299) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(1.6039675) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-2.0310903) q[3];
sx q[3];
rz(2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9017225) q[0];
sx q[0];
rz(-0.19342315) q[0];
sx q[0];
rz(1.1892009) q[0];
rz(0.72422782) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(1.6667746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14872197) q[0];
sx q[0];
rz(-1.9970511) q[0];
sx q[0];
rz(2.5919586) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15691136) q[2];
sx q[2];
rz(-2.3381066) q[2];
sx q[2];
rz(0.29565865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9660419) q[1];
sx q[1];
rz(-1.2331687) q[1];
sx q[1];
rz(-2.8199944) q[1];
rz(-pi) q[2];
rz(1.0737562) q[3];
sx q[3];
rz(-1.1599564) q[3];
sx q[3];
rz(-0.88676329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38179794) q[2];
sx q[2];
rz(-2.2319824) q[2];
sx q[2];
rz(-2.6749715) q[2];
rz(-1.7032547) q[3];
sx q[3];
rz(-1.8432901) q[3];
sx q[3];
rz(-2.2013825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.71317116) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(0.010350479) q[0];
rz(-1.6185282) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(-1.3759605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3165163) q[0];
sx q[0];
rz(-0.42046558) q[0];
sx q[0];
rz(-2.1964873) q[0];
x q[1];
rz(2.9508123) q[2];
sx q[2];
rz(-1.6990802) q[2];
sx q[2];
rz(0.14036638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37782323) q[1];
sx q[1];
rz(-1.0727912) q[1];
sx q[1];
rz(1.3643907) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67988396) q[3];
sx q[3];
rz(-2.2671642) q[3];
sx q[3];
rz(-1.0395466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4202262) q[2];
sx q[2];
rz(-2.0927636) q[2];
sx q[2];
rz(-2.4143207) q[2];
rz(0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0338243) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(-0.4069826) q[0];
rz(0.96626967) q[1];
sx q[1];
rz(-2.2844908) q[1];
sx q[1];
rz(0.13872096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1099595) q[0];
sx q[0];
rz(-3.0199625) q[0];
sx q[0];
rz(-1.2741791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66731668) q[2];
sx q[2];
rz(-0.34798589) q[2];
sx q[2];
rz(-1.1225357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87200621) q[1];
sx q[1];
rz(-1.2254189) q[1];
sx q[1];
rz(-0.38265444) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9646264) q[3];
sx q[3];
rz(-0.37016777) q[3];
sx q[3];
rz(0.98370508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8653284) q[2];
sx q[2];
rz(-1.7636969) q[2];
sx q[2];
rz(0.31065568) q[2];
rz(1.8336512) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(0.37180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0887611) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(-0.088223591) q[0];
rz(3.131033) q[1];
sx q[1];
rz(-1.5225531) q[1];
sx q[1];
rz(-2.6710076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6664526) q[0];
sx q[0];
rz(-1.0658403) q[0];
sx q[0];
rz(0.83374896) q[0];
rz(-pi) q[1];
rz(-1.3731277) q[2];
sx q[2];
rz(-0.77357793) q[2];
sx q[2];
rz(2.8171223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7631665) q[1];
sx q[1];
rz(-1.6554804) q[1];
sx q[1];
rz(1.3099109) q[1];
rz(-1.9235189) q[3];
sx q[3];
rz(-1.439487) q[3];
sx q[3];
rz(-2.2985947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(1.0799705) q[2];
rz(0.88932577) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(-1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.476986) q[0];
sx q[0];
rz(-0.38170686) q[0];
sx q[0];
rz(-2.3241924) q[0];
rz(-2.3951702) q[1];
sx q[1];
rz(-0.67459977) q[1];
sx q[1];
rz(0.93961632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25519263) q[0];
sx q[0];
rz(-1.5536683) q[0];
sx q[0];
rz(3.1411132) q[0];
x q[1];
rz(-1.290326) q[2];
sx q[2];
rz(-1.2277516) q[2];
sx q[2];
rz(-1.5743992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92781237) q[1];
sx q[1];
rz(-2.8646788) q[1];
sx q[1];
rz(2.0109948) q[1];
rz(1.1577179) q[3];
sx q[3];
rz(-1.7677757) q[3];
sx q[3];
rz(-2.911681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8115936) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(0.3453556) q[2];
rz(1.5251478) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(-0.42720544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.78373194) q[0];
sx q[0];
rz(-1.9851728) q[0];
sx q[0];
rz(0.082948908) q[0];
rz(2.9792765) q[1];
sx q[1];
rz(-1.6971308) q[1];
sx q[1];
rz(0.69581318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37980027) q[0];
sx q[0];
rz(-1.6125935) q[0];
sx q[0];
rz(0.27970741) q[0];
x q[1];
rz(-3.1386069) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(0.35831991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2999794) q[1];
sx q[1];
rz(-1.1265973) q[1];
sx q[1];
rz(1.5496553) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8105177) q[3];
sx q[3];
rz(-0.73270117) q[3];
sx q[3];
rz(0.98950451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7210228) q[2];
sx q[2];
rz(-0.15779725) q[2];
sx q[2];
rz(-1.921462) q[2];
rz(-2.5681791) q[3];
sx q[3];
rz(-1.3308595) q[3];
sx q[3];
rz(-2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.2777916) q[0];
sx q[0];
rz(-1.5024804) q[0];
sx q[0];
rz(-0.13308751) q[0];
rz(-3.1387023) q[1];
sx q[1];
rz(-0.4160226) q[1];
sx q[1];
rz(-0.70152534) q[1];
rz(1.3221424) q[2];
sx q[2];
rz(-1.7380309) q[2];
sx q[2];
rz(-2.0792014) q[2];
rz(3.0843432) q[3];
sx q[3];
rz(-2.3404239) q[3];
sx q[3];
rz(-0.41062582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
