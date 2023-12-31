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
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85672985) q[0];
sx q[0];
rz(-2.3861109) q[0];
sx q[0];
rz(1.5289686) q[0];
x q[1];
rz(0.062866048) q[2];
sx q[2];
rz(-2.2241484) q[2];
sx q[2];
rz(-1.9444998) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4418728) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-2.2295879) q[1];
x q[2];
rz(-0.05539031) q[3];
sx q[3];
rz(-0.4007383) q[3];
sx q[3];
rz(-1.8000359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(-0.43757004) q[0];
rz(-0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(0.22110573) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9843922) q[0];
sx q[0];
rz(-1.385014) q[0];
sx q[0];
rz(-1.1093344) q[0];
rz(-pi) q[1];
rz(1.8684623) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(0.01089451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088614956) q[1];
sx q[1];
rz(-2.6664475) q[1];
sx q[1];
rz(0.26762025) q[1];
rz(2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-0.75769889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706447) q[0];
sx q[0];
rz(-0.62951127) q[0];
sx q[0];
rz(-2.5802617) q[0];
rz(-pi) q[1];
rz(-0.65501113) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(0.26817817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.12049) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(2.2282269) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4593616) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.3285332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.3550718) q[0];
rz(-pi) q[1];
rz(-2.0650234) q[2];
sx q[2];
rz(-2.2708714) q[2];
sx q[2];
rz(1.9698524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7698344) q[1];
sx q[1];
rz(-1.8375988) q[1];
sx q[1];
rz(0.00053243551) q[1];
rz(-pi) q[2];
rz(-3.0321211) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(-2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(-1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(0.50450528) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.757526) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(1.9996044) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66820504) q[2];
sx q[2];
rz(-0.56958157) q[2];
sx q[2];
rz(2.6936206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0658768) q[1];
sx q[1];
rz(-1.2677691) q[1];
sx q[1];
rz(0.24995835) q[1];
x q[2];
rz(2.2284331) q[3];
sx q[3];
rz(-2.1214161) q[3];
sx q[3];
rz(0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(-1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
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
rz(-pi) q[1];
x q[1];
rz(1.2994453) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(-0.041610418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.530045) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(2.8987315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8998428) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-0.25776699) q[0];
sx q[0];
rz(-1.2447312) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6799913) q[2];
sx q[2];
rz(-1.9399376) q[2];
sx q[2];
rz(2.387407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(-1.741239) q[1];
rz(-0.93382436) q[3];
sx q[3];
rz(-1.0855506) q[3];
sx q[3];
rz(-2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-2.846472) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6587257) q[0];
sx q[0];
rz(-1.3728766) q[0];
sx q[0];
rz(-2.0107962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60111945) q[2];
sx q[2];
rz(-1.6496611) q[2];
sx q[2];
rz(2.2476946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1845491) q[1];
sx q[1];
rz(-2.7126985) q[1];
sx q[1];
rz(-1.2317608) q[1];
rz(-2.6929018) q[3];
sx q[3];
rz(-2.7033227) q[3];
sx q[3];
rz(-1.9697619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593489) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2851965) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(-1.517209) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8224045) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(-0.14949456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1907562) q[1];
sx q[1];
rz(-1.3637929) q[1];
sx q[1];
rz(2.2224269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82448126) q[3];
sx q[3];
rz(-1.609625) q[3];
sx q[3];
rz(2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(-2.8163731) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(2.7744055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(0.048006417) q[0];
rz(1.2904097) q[2];
sx q[2];
rz(-0.9624316) q[2];
sx q[2];
rz(0.3384564) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2328408) q[1];
sx q[1];
rz(-2.4526261) q[1];
sx q[1];
rz(-2.8244551) q[1];
rz(-0.92948593) q[3];
sx q[3];
rz(-1.9668285) q[3];
sx q[3];
rz(-3.0505153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29522482) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(2.8425343) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(0.085061442) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
