OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(-2.535948) q[1];
sx q[1];
rz(-0.65094596) q[1];
sx q[1];
rz(0.63408607) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91416042) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(3.1022275) q[0];
x q[1];
rz(-1.4889105) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(1.8413078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69971985) q[1];
sx q[1];
rz(-2.3708214) q[1];
sx q[1];
rz(-0.91200478) q[1];
rz(1.5942469) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32199931) q[0];
sx q[0];
rz(-2.0237192) q[0];
sx q[0];
rz(-0.20690147) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2220076) q[2];
sx q[2];
rz(-1.6752599) q[2];
sx q[2];
rz(-1.2806569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21048966) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(-1.4355852) q[1];
rz(-0.50237327) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(-2.3838938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0962778) q[0];
sx q[0];
rz(-2.0925539) q[0];
sx q[0];
rz(-1.2009215) q[0];
x q[1];
rz(2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(-0.26817817) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1670907) q[1];
sx q[1];
rz(-0.78364148) q[1];
sx q[1];
rz(-2.2553315) q[1];
rz(-pi) q[2];
rz(1.682231) q[3];
sx q[3];
rz(-0.93524869) q[3];
sx q[3];
rz(-3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.8130594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7610361) q[0];
sx q[0];
rz(-1.7839) q[0];
sx q[0];
rz(-2.9831432) q[0];
x q[1];
rz(0.76339108) q[2];
sx q[2];
rz(-1.9420468) q[2];
sx q[2];
rz(0.064918092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36973876) q[1];
sx q[1];
rz(-2.8747897) q[1];
sx q[1];
rz(-1.5727444) q[1];
rz(-1.7926746) q[3];
sx q[3];
rz(-1.4639877) q[3];
sx q[3];
rz(2.2282003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(-1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-2.6370874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.757526) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(1.1419883) q[0];
x q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(1.2010241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0658768) q[1];
sx q[1];
rz(-1.2677691) q[1];
sx q[1];
rz(2.8916343) q[1];
x q[2];
rz(-0.78318627) q[3];
sx q[3];
rz(-2.3110356) q[3];
sx q[3];
rz(-2.1577948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(-1.6882287) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(2.0507623) q[0];
rz(-0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-2.9528023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8180346) q[0];
sx q[0];
rz(-1.3209045) q[0];
sx q[0];
rz(2.7701993) q[0];
rz(-pi) q[1];
rz(1.2994453) q[2];
sx q[2];
rz(-1.4704629) q[2];
sx q[2];
rz(-3.0999822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.530045) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(0.24286119) q[1];
rz(-pi) q[2];
rz(-1.7536229) q[3];
sx q[3];
rz(-1.8086686) q[3];
sx q[3];
rz(2.6423303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.7623998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3430816) q[0];
sx q[0];
rz(-1.3268952) q[0];
sx q[0];
rz(-3.0573465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1629282) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-2.5025764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0849689) q[1];
sx q[1];
rz(-1.2722204) q[1];
sx q[1];
rz(1.741239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2961388) q[3];
sx q[3];
rz(-0.77973706) q[3];
sx q[3];
rz(-2.9065135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.032701187) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(-0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(0.76422894) q[0];
rz(0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.2932628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9613567) q[0];
sx q[0];
rz(-1.1399674) q[0];
sx q[0];
rz(2.9234617) q[0];
rz(-pi) q[1];
rz(0.60111945) q[2];
sx q[2];
rz(-1.4919315) q[2];
sx q[2];
rz(0.89389801) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.075899374) q[1];
sx q[1];
rz(-1.432044) q[1];
sx q[1];
rz(1.1636415) q[1];
rz(1.3702277) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(-1.6605103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-2.5602706) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2851965) q[0];
sx q[0];
rz(-2.6515793) q[0];
sx q[0];
rz(1.6243837) q[0];
x q[1];
rz(-1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(-0.21208866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-0.25823621) q[1];
x q[2];
rz(-3.0887293) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(2.7744055) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31736483) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-1.2904097) q[2];
sx q[2];
rz(-0.9624316) q[2];
sx q[2];
rz(-0.3384564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2314184) q[1];
sx q[1];
rz(-1.3712198) q[1];
sx q[1];
rz(-2.477596) q[1];
rz(-0.96079798) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-1.4916186) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(-2.9885837) q[3];
sx q[3];
rz(-0.85481337) q[3];
sx q[3];
rz(-0.92683642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
