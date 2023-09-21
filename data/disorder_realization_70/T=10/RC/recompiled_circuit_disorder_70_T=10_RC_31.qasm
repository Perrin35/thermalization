OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(3.7925386) q[1];
sx q[1];
rz(8.7906919) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91416042) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(-3.1022275) q[0];
x q[1];
rz(3.0787266) q[2];
sx q[2];
rz(-0.91744423) q[2];
sx q[2];
rz(-1.9444998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0178535) q[1];
sx q[1];
rz(-2.1542319) q[1];
sx q[1];
rz(2.6052193) q[1];
x q[2];
rz(-1.5942469) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.401702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32199931) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(0.20690147) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11110335) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(-2.8135516) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.088614956) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(0.26762025) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17640555) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(-2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0962778) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.2009215) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(0.26817817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97450199) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(-2.2553315) q[1];
rz(-0.14962872) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(3.0444037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(-1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(1.8130594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.3550718) q[0];
rz(-2.0650234) q[2];
sx q[2];
rz(-0.87072125) q[2];
sx q[2];
rz(1.1717403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1991785) q[1];
sx q[1];
rz(-1.5713099) q[1];
sx q[1];
rz(-1.3039939) q[1];
x q[2];
rz(3.0321211) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(1.1070586) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5915247) q[0];
sx q[0];
rz(-1.3378157) q[0];
sx q[0];
rz(1.0250807) q[0];
rz(-pi) q[1];
rz(0.66820504) q[2];
sx q[2];
rz(-0.56958157) q[2];
sx q[2];
rz(-0.44797209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(2.2401287) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4818146) q[3];
sx q[3];
rz(-1.0228844) q[3];
sx q[3];
rz(-1.9632615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.323558) q[0];
sx q[0];
rz(-1.3209045) q[0];
sx q[0];
rz(0.37139335) q[0];
x q[1];
rz(1.8421474) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(3.0999822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97555679) q[1];
sx q[1];
rz(-1.8131201) q[1];
sx q[1];
rz(1.6386599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24174989) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
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
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(1.8968614) q[0];
rz(2.6799913) q[2];
sx q[2];
rz(-1.9399376) q[2];
sx q[2];
rz(-0.75418562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52723253) q[1];
sx q[1];
rz(-2.7990606) q[1];
sx q[1];
rz(-0.50369461) q[1];
x q[2];
rz(-0.84545387) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11809764) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(2.3773637) q[0];
rz(3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30771502) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(1.1307554) q[0];
rz(-2.5404732) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(-0.89389801) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5871208) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(0.15092571) q[1];
rz(2.6929018) q[3];
sx q[3];
rz(-2.7033227) q[3];
sx q[3];
rz(-1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.8539799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33289136) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(1.0813792) q[0];
rz(1.7436696) q[2];
sx q[2];
rz(-1.0908974) q[2];
sx q[2];
rz(0.21208866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95083642) q[1];
sx q[1];
rz(-1.3637929) q[1];
sx q[1];
rz(-2.2224269) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0887293) q[3];
sx q[3];
rz(-2.3164146) q[3];
sx q[3];
rz(0.70096522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.256475) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(0.36718711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31736483) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-pi) q[1];
rz(-0.37819241) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(-2.3364002) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6348833) q[1];
sx q[1];
rz(-0.92223972) q[1];
sx q[1];
rz(1.3193921) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96079798) q[3];
sx q[3];
rz(-0.73878091) q[3];
sx q[3];
rz(-1.9566386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-1.6499741) q[2];
sx q[2];
rz(-0.82293246) q[2];
sx q[2];
rz(-0.056597829) q[2];
rz(-1.3973665) q[3];
sx q[3];
rz(-2.4122824) q[3];
sx q[3];
rz(-1.1576049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];