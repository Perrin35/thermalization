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
rz(0.68501002) q[0];
sx q[0];
rz(-2.4497439) q[0];
sx q[0];
rz(0.43842167) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(-2.9000403) q[1];
sx q[1];
rz(-0.99689364) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55132574) q[0];
sx q[0];
rz(-2.9661313) q[0];
sx q[0];
rz(2.0911123) q[0];
rz(0.27049944) q[2];
sx q[2];
rz(-1.6093614) q[2];
sx q[2];
rz(-0.46854101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6899851) q[1];
sx q[1];
rz(-1.8120736) q[1];
sx q[1];
rz(-1.1651768) q[1];
rz(-0.23714692) q[3];
sx q[3];
rz(-1.3714694) q[3];
sx q[3];
rz(-1.7502427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3991656) q[2];
sx q[2];
rz(-2.4319067) q[2];
sx q[2];
rz(-0.7134552) q[2];
rz(2.3136638) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(0.92638612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543318) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(1.2607505) q[0];
rz(0.24185355) q[1];
sx q[1];
rz(-1.8459903) q[1];
sx q[1];
rz(-0.58580011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59585786) q[0];
sx q[0];
rz(-2.7867466) q[0];
sx q[0];
rz(-1.097358) q[0];
rz(-0.42795534) q[2];
sx q[2];
rz(-1.5628067) q[2];
sx q[2];
rz(-2.7939579) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7603078) q[1];
sx q[1];
rz(-0.80712748) q[1];
sx q[1];
rz(1.155608) q[1];
x q[2];
rz(1.582566) q[3];
sx q[3];
rz(-1.571072) q[3];
sx q[3];
rz(-2.2451412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62023097) q[2];
sx q[2];
rz(-2.5197881) q[2];
sx q[2];
rz(0.28548959) q[2];
rz(2.1776958) q[3];
sx q[3];
rz(-0.6670835) q[3];
sx q[3];
rz(0.17531659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6730839) q[0];
sx q[0];
rz(-0.54247576) q[0];
sx q[0];
rz(2.8520404) q[0];
rz(0.30329224) q[1];
sx q[1];
rz(-1.8692317) q[1];
sx q[1];
rz(1.2264651) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2305812) q[0];
sx q[0];
rz(-0.93114656) q[0];
sx q[0];
rz(-0.53627642) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3704088) q[2];
sx q[2];
rz(-2.5287712) q[2];
sx q[2];
rz(-1.0072264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1947866) q[1];
sx q[1];
rz(-2.1420519) q[1];
sx q[1];
rz(2.8033546) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3555894) q[3];
sx q[3];
rz(-1.6542098) q[3];
sx q[3];
rz(2.047603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44846416) q[2];
sx q[2];
rz(-0.99702865) q[2];
sx q[2];
rz(0.86359751) q[2];
rz(0.10330769) q[3];
sx q[3];
rz(-1.3477252) q[3];
sx q[3];
rz(3.0062655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0057664) q[0];
sx q[0];
rz(-2.5209881) q[0];
sx q[0];
rz(-0.045106877) q[0];
rz(2.3182484) q[1];
sx q[1];
rz(-2.7594559) q[1];
sx q[1];
rz(-2.8270922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5161834) q[0];
sx q[0];
rz(-0.66358951) q[0];
sx q[0];
rz(-0.13701464) q[0];
rz(-2.7240924) q[2];
sx q[2];
rz(-1.0733777) q[2];
sx q[2];
rz(0.50499798) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.040995251) q[1];
sx q[1];
rz(-2.0597337) q[1];
sx q[1];
rz(1.3616427) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7847332) q[3];
sx q[3];
rz(-1.2535411) q[3];
sx q[3];
rz(-2.2283613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0752252) q[2];
sx q[2];
rz(-1.4868569) q[2];
sx q[2];
rz(-2.5328947) q[2];
rz(1.4501976) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(2.6628185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64959127) q[0];
sx q[0];
rz(-0.86968017) q[0];
sx q[0];
rz(-3.0539404) q[0];
rz(1.2610669) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(2.2671949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6325272) q[0];
sx q[0];
rz(-1.5380926) q[0];
sx q[0];
rz(2.2421809) q[0];
rz(-pi) q[1];
rz(2.1615209) q[2];
sx q[2];
rz(-1.5993759) q[2];
sx q[2];
rz(-2.9977552) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4147784) q[1];
sx q[1];
rz(-2.2764479) q[1];
sx q[1];
rz(-2.3470641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0239065) q[3];
sx q[3];
rz(-1.3804169) q[3];
sx q[3];
rz(3.0741907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38292357) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(0.60302889) q[2];
rz(-0.99271071) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(0.64811289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25931609) q[0];
sx q[0];
rz(-0.74463212) q[0];
sx q[0];
rz(-0.69389206) q[0];
rz(-0.71677417) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(-2.1778291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0345273) q[0];
sx q[0];
rz(-2.0052768) q[0];
sx q[0];
rz(-1.700842) q[0];
rz(-pi) q[1];
rz(-2.9593857) q[2];
sx q[2];
rz(-1.7413531) q[2];
sx q[2];
rz(-2.272416) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.239526) q[1];
sx q[1];
rz(-2.0888302) q[1];
sx q[1];
rz(0.43398989) q[1];
x q[2];
rz(3.0984466) q[3];
sx q[3];
rz(-0.57278297) q[3];
sx q[3];
rz(-1.6067895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-2.852735) q[2];
sx q[2];
rz(3.0333983) q[2];
rz(-0.9555971) q[3];
sx q[3];
rz(-0.038766131) q[3];
sx q[3];
rz(-0.55967104) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74681246) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(0.10699233) q[0];
rz(-0.50210285) q[1];
sx q[1];
rz(-1.6306449) q[1];
sx q[1];
rz(-0.14920251) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82115) q[0];
sx q[0];
rz(-1.4198729) q[0];
sx q[0];
rz(1.3942777) q[0];
x q[1];
rz(0.72609857) q[2];
sx q[2];
rz(-0.70961414) q[2];
sx q[2];
rz(0.47898705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1625357) q[1];
sx q[1];
rz(-1.874089) q[1];
sx q[1];
rz(2.3065662) q[1];
rz(-1.427569) q[3];
sx q[3];
rz(-1.9960072) q[3];
sx q[3];
rz(0.79637209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8788098) q[2];
sx q[2];
rz(-2.171319) q[2];
sx q[2];
rz(0.61857569) q[2];
rz(0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(-0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743643) q[0];
sx q[0];
rz(-1.3716797) q[0];
sx q[0];
rz(0.57269639) q[0];
rz(-2.122208) q[1];
sx q[1];
rz(-2.9044594) q[1];
sx q[1];
rz(-3.0651029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50216) q[0];
sx q[0];
rz(-2.0973839) q[0];
sx q[0];
rz(2.9798085) q[0];
rz(-pi) q[1];
rz(1.3504845) q[2];
sx q[2];
rz(-1.8615842) q[2];
sx q[2];
rz(2.5875768) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98497154) q[1];
sx q[1];
rz(-0.67423361) q[1];
sx q[1];
rz(-1.3944425) q[1];
x q[2];
rz(-0.88249607) q[3];
sx q[3];
rz(-1.9089437) q[3];
sx q[3];
rz(-0.32143092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1028334) q[2];
sx q[2];
rz(-0.88492727) q[2];
sx q[2];
rz(-2.6494675) q[2];
rz(-0.37880185) q[3];
sx q[3];
rz(-0.4327966) q[3];
sx q[3];
rz(0.88703275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4019796) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(0.42994764) q[0];
rz(-2.6509189) q[1];
sx q[1];
rz(-0.48777598) q[1];
sx q[1];
rz(1.0027764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802964) q[0];
sx q[0];
rz(-1.8334098) q[0];
sx q[0];
rz(-2.0735334) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.241469) q[2];
sx q[2];
rz(-1.5104093) q[2];
sx q[2];
rz(-2.4886069) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44081894) q[1];
sx q[1];
rz(-0.46844581) q[1];
sx q[1];
rz(3.0354795) q[1];
x q[2];
rz(-1.1100889) q[3];
sx q[3];
rz(-2.2177601) q[3];
sx q[3];
rz(-1.4972735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3873202) q[2];
sx q[2];
rz(-0.60101271) q[2];
sx q[2];
rz(2.4499272) q[2];
rz(0.12868853) q[3];
sx q[3];
rz(-1.5494989) q[3];
sx q[3];
rz(-2.5531829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9733031) q[0];
sx q[0];
rz(-3.0423218) q[0];
sx q[0];
rz(-0.6231935) q[0];
rz(-0.19459952) q[1];
sx q[1];
rz(-1.1293026) q[1];
sx q[1];
rz(-0.87619877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3850383) q[0];
sx q[0];
rz(-1.6858188) q[0];
sx q[0];
rz(-0.043941078) q[0];
x q[1];
rz(-0.16531971) q[2];
sx q[2];
rz(-2.39318) q[2];
sx q[2];
rz(0.95647631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9651803) q[1];
sx q[1];
rz(-1.1438055) q[1];
sx q[1];
rz(1.4216485) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.629452) q[3];
sx q[3];
rz(-0.63203963) q[3];
sx q[3];
rz(-1.2495981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96560043) q[2];
sx q[2];
rz(-0.24933641) q[2];
sx q[2];
rz(-0.636379) q[2];
rz(-1.6217568) q[3];
sx q[3];
rz(-0.88080019) q[3];
sx q[3];
rz(0.22708587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31430055) q[0];
sx q[0];
rz(-1.5999595) q[0];
sx q[0];
rz(-0.40524361) q[0];
rz(1.7397407) q[1];
sx q[1];
rz(-1.6442465) q[1];
sx q[1];
rz(0.13784611) q[1];
rz(-2.9707303) q[2];
sx q[2];
rz(-2.5124585) q[2];
sx q[2];
rz(-0.72825904) q[2];
rz(-0.53608175) q[3];
sx q[3];
rz(-1.5403219) q[3];
sx q[3];
rz(2.2684682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
