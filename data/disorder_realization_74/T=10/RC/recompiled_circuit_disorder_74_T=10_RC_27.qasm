OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9292944) q[0];
sx q[0];
rz(-1.1120783) q[0];
sx q[0];
rz(-0.32231583) q[0];
rz(-pi) q[1];
rz(2.870954) q[2];
sx q[2];
rz(-0.025279609) q[2];
sx q[2];
rz(1.4852922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9413486) q[1];
sx q[1];
rz(-1.3625047) q[1];
sx q[1];
rz(-1.3079446) q[1];
x q[2];
rz(-0.07006499) q[3];
sx q[3];
rz(-1.8383887) q[3];
sx q[3];
rz(-0.90581264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(-0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-2.8285817) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(1.82812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.388893) q[0];
sx q[0];
rz(-1.3503195) q[0];
sx q[0];
rz(-1.3541143) q[0];
rz(-pi) q[1];
rz(0.30828373) q[2];
sx q[2];
rz(-2.0364025) q[2];
sx q[2];
rz(-0.75454933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7105512) q[1];
sx q[1];
rz(-1.6574873) q[1];
sx q[1];
rz(-2.980742) q[1];
rz(-3.0659862) q[3];
sx q[3];
rz(-1.6763655) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(-2.4296956) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-0.16608873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8649193) q[0];
sx q[0];
rz(-1.2827946) q[0];
sx q[0];
rz(-2.9898781) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(-2.6454676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3661763) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(2.612545) q[1];
rz(2.283964) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(1.5989725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(3.0825901) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(0.44149533) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(-0.77484432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7168717) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(-1.5543235) q[0];
x q[1];
rz(1.5261126) q[2];
sx q[2];
rz(-1.9384906) q[2];
sx q[2];
rz(-1.1553264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1244509) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(3.0567398) q[1];
x q[2];
rz(0.99153783) q[3];
sx q[3];
rz(-2.8354128) q[3];
sx q[3];
rz(1.7508208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775979) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(-0.60896215) q[1];
sx q[1];
rz(-2.6700171) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5126702) q[0];
sx q[0];
rz(-1.218601) q[0];
sx q[0];
rz(0.098350071) q[0];
x q[1];
rz(-0.92685076) q[2];
sx q[2];
rz(-1.1139718) q[2];
sx q[2];
rz(-0.18075519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4957461) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(0.87162019) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34541901) q[3];
sx q[3];
rz(-1.0343699) q[3];
sx q[3];
rz(-0.1766583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.057377664) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058379563) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-1.0169792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4646343) q[0];
sx q[0];
rz(-2.5559253) q[0];
sx q[0];
rz(2.2023489) q[0];
rz(-pi) q[1];
rz(0.20118841) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(-2.3692998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79343866) q[1];
sx q[1];
rz(-0.70421709) q[1];
sx q[1];
rz(-2.3401295) q[1];
x q[2];
rz(1.059504) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0680925) q[0];
sx q[0];
rz(-1.4788027) q[0];
sx q[0];
rz(1.5324701) q[0];
rz(0.90791038) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(-2.8223035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19983229) q[1];
sx q[1];
rz(-0.70711771) q[1];
sx q[1];
rz(1.1862399) q[1];
x q[2];
rz(0.089902417) q[3];
sx q[3];
rz(-1.1635167) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17657875) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(-2.2289508) q[0];
rz(1.2572844) q[2];
sx q[2];
rz(-1.0672617) q[2];
sx q[2];
rz(1.7538278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9600535) q[1];
sx q[1];
rz(-1.0078733) q[1];
sx q[1];
rz(1.3905418) q[1];
rz(2.7301634) q[3];
sx q[3];
rz(-2.5860388) q[3];
sx q[3];
rz(3.0761449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2239969) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.6329637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471837) q[0];
sx q[0];
rz(-1.7315454) q[0];
sx q[0];
rz(0.87662351) q[0];
x q[1];
rz(1.6998859) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(2.8512851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3907527) q[1];
sx q[1];
rz(-0.62493338) q[1];
sx q[1];
rz(1.0052488) q[1];
rz(-1.6375332) q[3];
sx q[3];
rz(-0.70422322) q[3];
sx q[3];
rz(-0.080554068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.7601097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.5230595) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6818468) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(-1.2088838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6227464) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-2.3153789) q[1];
rz(-1.8990717) q[3];
sx q[3];
rz(-0.74088135) q[3];
sx q[3];
rz(1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(3.1141282) q[3];
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
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(1.7455208) q[2];
sx q[2];
rz(-2.8833485) q[2];
sx q[2];
rz(-1.9042518) q[2];
rz(-0.14470312) q[3];
sx q[3];
rz(-2.6270513) q[3];
sx q[3];
rz(-2.4270121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
