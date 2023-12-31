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
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68361002) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(0.81575127) q[0];
x q[1];
rz(1.4889105) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(1.3002849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4418728) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-0.91200478) q[1];
rz(-pi) q[2];
rz(-1.5473458) q[3];
sx q[3];
rz(-1.970885) q[3];
sx q[3];
rz(1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32199931) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(-2.9346912) q[0];
rz(-3.0304893) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(0.32804104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.088614956) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(-0.26762025) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17640555) q[3];
sx q[3];
rz(-2.6325912) q[3];
sx q[3];
rz(-2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(2.5533) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(2.3838938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43514566) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
x q[1];
rz(2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(-0.26817817) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0231087) q[1];
sx q[1];
rz(-0.99220905) q[1];
sx q[1];
rz(2.5793377) q[1];
rz(-pi) q[2];
x q[2];
rz(1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(0.20382717) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(1.3285332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(-1.7865208) q[0];
rz(-pi) q[1];
rz(2.3782016) q[2];
sx q[2];
rz(-1.9420468) q[2];
sx q[2];
rz(3.0766746) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7698344) q[1];
sx q[1];
rz(-1.3039939) q[1];
sx q[1];
rz(-0.00053243551) q[1];
rz(-1.1174326) q[3];
sx q[3];
rz(-0.24586596) q[3];
sx q[3];
rz(2.9256431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5915247) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(-1.0250807) q[0];
rz(-pi) q[1];
rz(0.66820504) q[2];
sx q[2];
rz(-2.5720111) q[2];
sx q[2];
rz(-2.6936206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5711172) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(1.8829324) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65977804) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(-1.9632615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8180346) q[0];
sx q[0];
rz(-1.3209045) q[0];
sx q[0];
rz(0.37139335) q[0];
rz(-pi) q[1];
rz(1.8421474) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(-0.041610418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61154762) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(-0.24286119) q[1];
rz(-pi) q[2];
rz(2.4981899) q[3];
sx q[3];
rz(-0.2989558) q[3];
sx q[3];
rz(-1.1645731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-1.0287372) q[2];
sx q[2];
rz(-1.8264654) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(1.3791929) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(1.8968614) q[0];
rz(0.71521476) q[2];
sx q[2];
rz(-2.5589802) q[2];
sx q[2];
rz(-1.4441393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56475793) q[1];
sx q[1];
rz(-1.4079637) q[1];
sx q[1];
rz(2.8388883) q[1];
rz(-pi) q[2];
rz(2.2961388) q[3];
sx q[3];
rz(-2.3618556) q[3];
sx q[3];
rz(-2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.8483298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30771502) q[0];
sx q[0];
rz(-0.47979646) q[0];
sx q[0];
rz(-2.0108372) q[0];
x q[1];
rz(-1.4752611) q[2];
sx q[2];
rz(-2.1697858) q[2];
sx q[2];
rz(2.4107188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5544719) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(2.9906669) q[1];
rz(-0.44869081) q[3];
sx q[3];
rz(-0.43827) q[3];
sx q[3];
rz(-1.9697619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(1.2506437) q[0];
rz(-2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8087013) q[0];
sx q[0];
rz(-1.5455855) q[0];
sx q[0];
rz(-1.0813792) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7436696) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(0.21208866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3567644) q[1];
sx q[1];
rz(-2.462466) q[1];
sx q[1];
rz(1.9041512) q[1];
rz(1.5136396) q[3];
sx q[3];
rz(-2.3944629) q[3];
sx q[3];
rz(-2.3627918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(-2.9750032) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-2.0064158) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(0.048006417) q[0];
x q[1];
rz(1.2904097) q[2];
sx q[2];
rz(-2.1791611) q[2];
sx q[2];
rz(-0.3384564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2314184) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(-2.477596) q[1];
x q[2];
rz(-2.2121067) q[3];
sx q[3];
rz(-1.9668285) q[3];
sx q[3];
rz(3.0505153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
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
rz(0.15300898) q[3];
sx q[3];
rz(-0.85481337) q[3];
sx q[3];
rz(-0.92683642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
