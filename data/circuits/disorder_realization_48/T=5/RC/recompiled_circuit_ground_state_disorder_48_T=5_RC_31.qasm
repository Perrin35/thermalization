OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(1.9999003) q[0];
sx q[0];
rz(6.5394149) q[0];
rz(-2.7645219) q[1];
sx q[1];
rz(-1.3426251) q[1];
sx q[1];
rz(-1.0599729) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621583) q[0];
sx q[0];
rz(-2.9720364) q[0];
sx q[0];
rz(0.79742278) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5358309) q[2];
sx q[2];
rz(-0.87858534) q[2];
sx q[2];
rz(-3.0193605) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6457451) q[1];
sx q[1];
rz(-1.331067) q[1];
sx q[1];
rz(-1.7868397) q[1];
x q[2];
rz(-1.6698631) q[3];
sx q[3];
rz(-1.7145673) q[3];
sx q[3];
rz(-3.037775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18067351) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.712435) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(3.070224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(-0.74547705) q[0];
rz(1.2019134) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(-2.1048996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1069477) q[0];
sx q[0];
rz(-1.4781524) q[0];
sx q[0];
rz(-1.9361467) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5971423) q[2];
sx q[2];
rz(-0.80412946) q[2];
sx q[2];
rz(-2.1955817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0526317) q[1];
sx q[1];
rz(-0.76238576) q[1];
sx q[1];
rz(-0.41684581) q[1];
rz(1.3535301) q[3];
sx q[3];
rz(-0.99743069) q[3];
sx q[3];
rz(2.2305302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1408954) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(-2.1689283) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.5001517) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(-2.3386173) q[0];
rz(-2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48939368) q[0];
sx q[0];
rz(-0.94871657) q[0];
sx q[0];
rz(-2.6724044) q[0];
rz(1.8016544) q[2];
sx q[2];
rz(-0.95536026) q[2];
sx q[2];
rz(-2.1417257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4068702) q[1];
sx q[1];
rz(-2.0592505) q[1];
sx q[1];
rz(3.0367756) q[1];
x q[2];
rz(2.8057116) q[3];
sx q[3];
rz(-1.349873) q[3];
sx q[3];
rz(2.9136052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.190072) q[2];
sx q[2];
rz(-2.0900487) q[2];
sx q[2];
rz(1.5104843) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-1.0534143) q[3];
sx q[3];
rz(-2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25201061) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(2.4776283) q[0];
rz(0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(1.0391611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1215724) q[0];
sx q[0];
rz(-1.5680321) q[0];
sx q[0];
rz(3.1317668) q[0];
rz(-pi) q[1];
rz(-1.9661994) q[2];
sx q[2];
rz(-1.3229473) q[2];
sx q[2];
rz(1.445766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0137456) q[1];
sx q[1];
rz(-1.3617853) q[1];
sx q[1];
rz(-1.0793988) q[1];
rz(-pi) q[2];
rz(0.19840513) q[3];
sx q[3];
rz(-2.6909749) q[3];
sx q[3];
rz(0.042447986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.7258464) q[2];
sx q[2];
rz(-3.0636129) q[2];
rz(1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(-1.0361205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2074821) q[0];
sx q[0];
rz(-2.9353862) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(1.9059937) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(1.3135757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.006047) q[0];
sx q[0];
rz(-2.1925547) q[0];
sx q[0];
rz(-2.9154587) q[0];
rz(-pi) q[1];
rz(-2.7817725) q[2];
sx q[2];
rz(-2.0096547) q[2];
sx q[2];
rz(0.0024589389) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2755679) q[1];
sx q[1];
rz(-2.0179664) q[1];
sx q[1];
rz(-0.42110301) q[1];
rz(-pi) q[2];
rz(0.4540654) q[3];
sx q[3];
rz(-2.3010074) q[3];
sx q[3];
rz(-0.24475748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95973394) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(-2.6971297) q[2];
rz(-1.2005165) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(-0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(-0.077202395) q[0];
rz(0.96241799) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(0.033800689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0173967) q[0];
sx q[0];
rz(-0.58308812) q[0];
sx q[0];
rz(-2.4840607) q[0];
rz(-pi) q[1];
rz(-2.5939213) q[2];
sx q[2];
rz(-1.593809) q[2];
sx q[2];
rz(-0.097824899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.046545854) q[1];
sx q[1];
rz(-2.425417) q[1];
sx q[1];
rz(-1.0821872) q[1];
x q[2];
rz(0.77671364) q[3];
sx q[3];
rz(-1.4956105) q[3];
sx q[3];
rz(-3.1326339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(-0.25924337) q[2];
rz(2.1974468) q[3];
sx q[3];
rz(-2.9278432) q[3];
sx q[3];
rz(1.8102144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0807226) q[0];
sx q[0];
rz(-2.935077) q[0];
sx q[0];
rz(0.61332214) q[0];
rz(-0.90336409) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(-0.14796743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4658303) q[0];
sx q[0];
rz(-1.0935385) q[0];
sx q[0];
rz(0.81264074) q[0];
rz(-pi) q[1];
x q[1];
rz(1.602722) q[2];
sx q[2];
rz(-1.1319926) q[2];
sx q[2];
rz(-1.1915468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92098713) q[1];
sx q[1];
rz(-1.1488394) q[1];
sx q[1];
rz(0.76442952) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5336125) q[3];
sx q[3];
rz(-2.5761009) q[3];
sx q[3];
rz(2.9786106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.938574) q[2];
sx q[2];
rz(-1.4116986) q[2];
sx q[2];
rz(2.1417248) q[2];
rz(-1.8262919) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(-0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(-1.8723764) q[1];
sx q[1];
rz(-2.3321584) q[1];
sx q[1];
rz(0.16407897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7456147) q[0];
sx q[0];
rz(-1.656129) q[0];
sx q[0];
rz(0.11767582) q[0];
x q[1];
rz(2.9226926) q[2];
sx q[2];
rz(-2.0563807) q[2];
sx q[2];
rz(-1.6644434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33262256) q[1];
sx q[1];
rz(-1.5676985) q[1];
sx q[1];
rz(-2.0634406) q[1];
rz(-2.9470575) q[3];
sx q[3];
rz(-0.92960581) q[3];
sx q[3];
rz(-1.7879144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8998731) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(-2.0797753) q[2];
rz(3.0070983) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(-3.0883446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4835085) q[0];
sx q[0];
rz(-2.4299419) q[0];
sx q[0];
rz(-2.9363976) q[0];
rz(-0.23458734) q[1];
sx q[1];
rz(-1.7154452) q[1];
sx q[1];
rz(0.1964143) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72442833) q[0];
sx q[0];
rz(-2.1311893) q[0];
sx q[0];
rz(-2.7743894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26650776) q[2];
sx q[2];
rz(-2.6015913) q[2];
sx q[2];
rz(-0.85734474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82746303) q[1];
sx q[1];
rz(-1.2478831) q[1];
sx q[1];
rz(2.9548378) q[1];
rz(-0.78077448) q[3];
sx q[3];
rz(-1.1765683) q[3];
sx q[3];
rz(0.0782644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3966169) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.1783925) q[2];
rz(0.34934238) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0013393764) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(0.70571357) q[0];
rz(-0.69752518) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(-0.90550214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4587935) q[0];
sx q[0];
rz(-1.1371277) q[0];
sx q[0];
rz(-1.5894029) q[0];
x q[1];
rz(2.7780813) q[2];
sx q[2];
rz(-1.8773407) q[2];
sx q[2];
rz(2.6078122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9230629) q[1];
sx q[1];
rz(-2.8273812) q[1];
sx q[1];
rz(0.15337069) q[1];
rz(1.9344744) q[3];
sx q[3];
rz(-1.26059) q[3];
sx q[3];
rz(-2.9077173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8991578) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(2.4230797) q[2];
rz(-0.68273035) q[3];
sx q[3];
rz(-1.9971137) q[3];
sx q[3];
rz(0.1040641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.318442) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(1.6406583) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(1.4529543) q[2];
sx q[2];
rz(-1.4284882) q[2];
sx q[2];
rz(1.535648) q[2];
rz(0.90632306) q[3];
sx q[3];
rz(-0.88574468) q[3];
sx q[3];
rz(1.1024324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
