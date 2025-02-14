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
rz(-1.4122352) q[0];
sx q[0];
rz(-0.5923624) q[0];
sx q[0];
rz(-1.2047729) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(-1.187721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43150768) q[0];
sx q[0];
rz(-1.0241226) q[0];
sx q[0];
rz(1.3923079) q[0];
rz(-0.6600136) q[2];
sx q[2];
rz(-2.4610991) q[2];
sx q[2];
rz(-1.6689491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7379308) q[1];
sx q[1];
rz(-1.523535) q[1];
sx q[1];
rz(-1.3980327) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0023469) q[3];
sx q[3];
rz(-1.2069824) q[3];
sx q[3];
rz(1.1195967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(1.0428693) q[2];
rz(0.045182191) q[3];
sx q[3];
rz(-0.16862814) q[3];
sx q[3];
rz(1.5359623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(-0.97852069) q[0];
rz(-1.9784031) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.8189583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90962553) q[0];
sx q[0];
rz(-2.2920906) q[0];
sx q[0];
rz(2.2366174) q[0];
rz(0.18206667) q[2];
sx q[2];
rz(-1.2918424) q[2];
sx q[2];
rz(-1.830991) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95056995) q[1];
sx q[1];
rz(-2.2132067) q[1];
sx q[1];
rz(0.2705785) q[1];
rz(-pi) q[2];
rz(0.17830542) q[3];
sx q[3];
rz(-2.4394991) q[3];
sx q[3];
rz(-2.3012379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3375552) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(-1.6328968) q[3];
sx q[3];
rz(-1.326606) q[3];
sx q[3];
rz(-1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22817336) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(-0.2359373) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.6940247) q[1];
sx q[1];
rz(2.6441914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.936393) q[0];
sx q[0];
rz(-2.8872364) q[0];
sx q[0];
rz(-0.78820552) q[0];
rz(-1.1549306) q[2];
sx q[2];
rz(-1.207104) q[2];
sx q[2];
rz(1.1689344) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.062089109) q[1];
sx q[1];
rz(-1.5928643) q[1];
sx q[1];
rz(-2.7479321) q[1];
x q[2];
rz(2.3373691) q[3];
sx q[3];
rz(-0.35103335) q[3];
sx q[3];
rz(1.1277175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(-2.2836397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-2.9486935) q[0];
rz(1.7861722) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(-0.17098175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4701914) q[0];
sx q[0];
rz(-2.6220052) q[0];
sx q[0];
rz(1.4524824) q[0];
x q[1];
rz(-2.4055482) q[2];
sx q[2];
rz(-1.9878329) q[2];
sx q[2];
rz(-1.2936178) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4857085) q[1];
sx q[1];
rz(-7*pi/15) q[1];
sx q[1];
rz(1.1397362) q[1];
rz(-pi) q[2];
rz(1.8754575) q[3];
sx q[3];
rz(-0.6864635) q[3];
sx q[3];
rz(-2.1368419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96044797) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-0.32285264) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(-2.3006181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5906931) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(2.2160231) q[0];
rz(0.12807056) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(0.099460348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8444237) q[0];
sx q[0];
rz(-2.4250829) q[0];
sx q[0];
rz(1.5708718) q[0];
x q[1];
rz(2.3771466) q[2];
sx q[2];
rz(-0.4808772) q[2];
sx q[2];
rz(2.1493727) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4872553) q[1];
sx q[1];
rz(-2.0981826) q[1];
sx q[1];
rz(-0.22589639) q[1];
rz(-pi) q[2];
rz(1.3252898) q[3];
sx q[3];
rz(-1.6485896) q[3];
sx q[3];
rz(0.93279749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(0.37364513) q[2];
rz(0.89961189) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(1.1348178) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(-2.5496971) q[0];
rz(-2.9369211) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(-0.051102292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44296023) q[0];
sx q[0];
rz(-1.763183) q[0];
sx q[0];
rz(1.7257742) q[0];
rz(-0.25801664) q[2];
sx q[2];
rz(-2.2126901) q[2];
sx q[2];
rz(0.24680617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7390427) q[1];
sx q[1];
rz(-1.879917) q[1];
sx q[1];
rz(1.9376631) q[1];
rz(-pi) q[2];
rz(-1.072364) q[3];
sx q[3];
rz(-1.0921156) q[3];
sx q[3];
rz(2.5414027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11703141) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(2.7545605) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36560202) q[0];
sx q[0];
rz(-0.87823534) q[0];
sx q[0];
rz(2.8940417) q[0];
rz(-1.4546825) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(-3.1051342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582583) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-2.3733632) q[0];
x q[1];
rz(-0.22561947) q[2];
sx q[2];
rz(-0.7157514) q[2];
sx q[2];
rz(-2.4168487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6853906) q[1];
sx q[1];
rz(-1.0674607) q[1];
sx q[1];
rz(-0.53629843) q[1];
x q[2];
rz(0.8730252) q[3];
sx q[3];
rz(-0.11618488) q[3];
sx q[3];
rz(-2.1538863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(-1.3521693) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9098814) q[0];
sx q[0];
rz(-0.24869643) q[0];
sx q[0];
rz(-1.6491718) q[0];
rz(0.17852783) q[1];
sx q[1];
rz(-1.8793722) q[1];
sx q[1];
rz(-0.55007225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290925) q[0];
sx q[0];
rz(-1.5508473) q[0];
sx q[0];
rz(1.602142) q[0];
x q[1];
rz(0.063385609) q[2];
sx q[2];
rz(-1.3235385) q[2];
sx q[2];
rz(-2.4154369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4143205) q[1];
sx q[1];
rz(-1.4106693) q[1];
sx q[1];
rz(-0.54267197) q[1];
x q[2];
rz(0.40591235) q[3];
sx q[3];
rz(-2.0024869) q[3];
sx q[3];
rz(1.622705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.4057691) q[2];
sx q[2];
rz(-2.3883635) q[2];
rz(-0.29427648) q[3];
sx q[3];
rz(-3.0615276) q[3];
sx q[3];
rz(-3.0729955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8681965) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(1.4519325) q[0];
rz(0.1952576) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(1.6498227) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378252) q[0];
sx q[0];
rz(-1.2653102) q[0];
sx q[0];
rz(2.6432493) q[0];
rz(-pi) q[1];
rz(-3.0259616) q[2];
sx q[2];
rz(-2.3732269) q[2];
sx q[2];
rz(2.7344784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10354708) q[1];
sx q[1];
rz(-2.9730151) q[1];
sx q[1];
rz(1.9650616) q[1];
rz(-2.9025495) q[3];
sx q[3];
rz(-0.028246917) q[3];
sx q[3];
rz(2.7735428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1008272) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(0.85774285) q[2];
rz(-1.6929251) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0070852) q[0];
sx q[0];
rz(-2.9869098) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(1.1231517) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(-0.060308594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6864844) q[0];
sx q[0];
rz(-1.9251072) q[0];
sx q[0];
rz(0.099781009) q[0];
x q[1];
rz(-3.0357828) q[2];
sx q[2];
rz(-2.1544837) q[2];
sx q[2];
rz(-0.90487827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8728208) q[1];
sx q[1];
rz(-2.4686681) q[1];
sx q[1];
rz(-3.0329143) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7785886) q[3];
sx q[3];
rz(-1.6771183) q[3];
sx q[3];
rz(0.3029938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-0.86468148) q[2];
sx q[2];
rz(-1.8314499) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.3430877) q[3];
sx q[3];
rz(-1.3069299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46351984) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(0.75640596) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(-1.5613772) q[2];
sx q[2];
rz(-1.1404372) q[2];
sx q[2];
rz(-2.5802317) q[2];
rz(-2.9587119) q[3];
sx q[3];
rz(-1.9149078) q[3];
sx q[3];
rz(-0.2855466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
