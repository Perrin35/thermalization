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
rz(0.19658495) q[0];
sx q[0];
rz(-2.6994446) q[0];
sx q[0];
rz(2.2801939) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(-0.035592508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1382768) q[0];
sx q[0];
rz(-0.39130032) q[0];
sx q[0];
rz(0.53404339) q[0];
x q[1];
rz(-2.8350949) q[2];
sx q[2];
rz(-1.674429) q[2];
sx q[2];
rz(0.95970861) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5000677) q[1];
sx q[1];
rz(-1.6065805) q[1];
sx q[1];
rz(-2.4149717) q[1];
rz(-0.017982131) q[3];
sx q[3];
rz(-2.2507022) q[3];
sx q[3];
rz(2.3462092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5888136) q[2];
sx q[2];
rz(-1.8401044) q[2];
sx q[2];
rz(-1.0944875) q[2];
rz(1.5158481) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(-0.14357963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0482408) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(0.37407237) q[0];
rz(2.2834942) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(1.0757974) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5257403) q[0];
sx q[0];
rz(-1.4676371) q[0];
sx q[0];
rz(0.10248689) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38041842) q[2];
sx q[2];
rz(-0.58092344) q[2];
sx q[2];
rz(-0.38512938) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.058671) q[1];
sx q[1];
rz(-1.6870912) q[1];
sx q[1];
rz(-2.7398159) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2298563) q[3];
sx q[3];
rz(-2.1359332) q[3];
sx q[3];
rz(-2.679058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51521236) q[2];
sx q[2];
rz(-2.7084646) q[2];
sx q[2];
rz(-2.058378) q[2];
rz(1.2927239) q[3];
sx q[3];
rz(-2.0079398) q[3];
sx q[3];
rz(-0.92866549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.885963) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-2.6166925) q[0];
rz(-1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(0.14579138) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12851579) q[0];
sx q[0];
rz(-1.6616964) q[0];
sx q[0];
rz(3.035782) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86051382) q[2];
sx q[2];
rz(-2.4053221) q[2];
sx q[2];
rz(-0.78303534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5151095) q[1];
sx q[1];
rz(-1.6609539) q[1];
sx q[1];
rz(-2.7075591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5691312) q[3];
sx q[3];
rz(-2.2160071) q[3];
sx q[3];
rz(2.0544586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24474457) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(-2.9761918) q[2];
rz(-2.2798955) q[3];
sx q[3];
rz(-0.69178897) q[3];
sx q[3];
rz(-0.19680463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2794613) q[0];
sx q[0];
rz(-2.617351) q[0];
sx q[0];
rz(-0.72823802) q[0];
rz(-0.32314745) q[1];
sx q[1];
rz(-2.1073982) q[1];
sx q[1];
rz(-2.1023777) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82074219) q[0];
sx q[0];
rz(-0.59950638) q[0];
sx q[0];
rz(2.3248424) q[0];
rz(0.93499281) q[2];
sx q[2];
rz(-1.7991009) q[2];
sx q[2];
rz(-1.993597) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8700712) q[1];
sx q[1];
rz(-1.0876552) q[1];
sx q[1];
rz(-1.6752233) q[1];
rz(-pi) q[2];
rz(2.3226854) q[3];
sx q[3];
rz(-2.120813) q[3];
sx q[3];
rz(0.14980042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6497961) q[2];
sx q[2];
rz(-2.0911262) q[2];
sx q[2];
rz(-2.6226079) q[2];
rz(-2.0670048) q[3];
sx q[3];
rz(-1.2858177) q[3];
sx q[3];
rz(-2.6463267) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062926) q[0];
sx q[0];
rz(-2.2195897) q[0];
sx q[0];
rz(-2.9700188) q[0];
rz(-1.0942787) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(-0.23460728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5253403) q[0];
sx q[0];
rz(-1.8484637) q[0];
sx q[0];
rz(3.0086631) q[0];
rz(0.50723664) q[2];
sx q[2];
rz(-1.6167069) q[2];
sx q[2];
rz(0.017233726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68069862) q[1];
sx q[1];
rz(-2.1354224) q[1];
sx q[1];
rz(-0.84565296) q[1];
rz(-1.9453641) q[3];
sx q[3];
rz(-2.2456333) q[3];
sx q[3];
rz(2.5684211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2285063) q[2];
sx q[2];
rz(-1.9501016) q[2];
sx q[2];
rz(1.7388434) q[2];
rz(-0.62471041) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(1.3431842) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45873555) q[0];
sx q[0];
rz(-3.1338437) q[0];
sx q[0];
rz(-0.32753456) q[0];
rz(3.1134743) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-0.99172529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1234489) q[0];
sx q[0];
rz(-0.94052343) q[0];
sx q[0];
rz(2.2176377) q[0];
x q[1];
rz(-0.4430094) q[2];
sx q[2];
rz(-2.8923419) q[2];
sx q[2];
rz(-0.64068782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8918482) q[1];
sx q[1];
rz(-1.5288903) q[1];
sx q[1];
rz(-0.81467198) q[1];
x q[2];
rz(3.087849) q[3];
sx q[3];
rz(-0.47131594) q[3];
sx q[3];
rz(1.1137247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.907054) q[2];
sx q[2];
rz(-1.3776642) q[2];
sx q[2];
rz(1.6592337) q[2];
rz(-1.0659069) q[3];
sx q[3];
rz(-1.9612471) q[3];
sx q[3];
rz(2.0445686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3858353) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(-0.29543153) q[0];
rz(-2.9040728) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(-0.74737731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4656671) q[0];
sx q[0];
rz(-2.3479778) q[0];
sx q[0];
rz(0.74489745) q[0];
rz(-pi) q[1];
rz(2.5146444) q[2];
sx q[2];
rz(-2.0184506) q[2];
sx q[2];
rz(0.7426151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44386266) q[1];
sx q[1];
rz(-1.2011114) q[1];
sx q[1];
rz(1.0401985) q[1];
rz(-0.57716863) q[3];
sx q[3];
rz(-0.62466972) q[3];
sx q[3];
rz(0.28055252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6682917) q[2];
sx q[2];
rz(-1.1857727) q[2];
sx q[2];
rz(2.8866923) q[2];
rz(0.21236803) q[3];
sx q[3];
rz(-0.86447159) q[3];
sx q[3];
rz(-0.00096850639) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85106987) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(-3.0176924) q[0];
rz(-2.2663785) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(-0.55878729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717455) q[0];
sx q[0];
rz(-2.5515243) q[0];
sx q[0];
rz(-1.1099932) q[0];
x q[1];
rz(-2.8577639) q[2];
sx q[2];
rz(-1.0887732) q[2];
sx q[2];
rz(-1.1555156) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81441036) q[1];
sx q[1];
rz(-2.5375527) q[1];
sx q[1];
rz(2.6534897) q[1];
x q[2];
rz(2.111192) q[3];
sx q[3];
rz(-0.70941389) q[3];
sx q[3];
rz(-1.7075001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.045804068) q[2];
sx q[2];
rz(-2.6498821) q[2];
sx q[2];
rz(0.72074214) q[2];
rz(-0.36302429) q[3];
sx q[3];
rz(-1.2237153) q[3];
sx q[3];
rz(1.629841) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1849798) q[0];
sx q[0];
rz(-2.9473801) q[0];
sx q[0];
rz(0.089056253) q[0];
rz(2.7748499) q[1];
sx q[1];
rz(-0.9042424) q[1];
sx q[1];
rz(1.6820224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37292591) q[0];
sx q[0];
rz(-0.67250508) q[0];
sx q[0];
rz(-1.1993394) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.681293) q[2];
sx q[2];
rz(-0.80901399) q[2];
sx q[2];
rz(-1.9495277) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9391528) q[1];
sx q[1];
rz(-1.9382361) q[1];
sx q[1];
rz(2.3177486) q[1];
rz(2.5404471) q[3];
sx q[3];
rz(-0.69128762) q[3];
sx q[3];
rz(2.9959681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2064994) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(-1.8915141) q[2];
rz(-1.2262723) q[3];
sx q[3];
rz(-2.5076702) q[3];
sx q[3];
rz(2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117689) q[0];
sx q[0];
rz(-1.1331695) q[0];
sx q[0];
rz(-1.1400219) q[0];
rz(2.5584768) q[1];
sx q[1];
rz(-1.6551599) q[1];
sx q[1];
rz(0.66782943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9662922) q[0];
sx q[0];
rz(-0.81476962) q[0];
sx q[0];
rz(0.45481429) q[0];
rz(-pi) q[1];
rz(-1.2905707) q[2];
sx q[2];
rz(-0.72151583) q[2];
sx q[2];
rz(1.0204329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17396233) q[1];
sx q[1];
rz(-0.47904992) q[1];
sx q[1];
rz(-1.1381989) q[1];
rz(2.0290501) q[3];
sx q[3];
rz(-0.83247165) q[3];
sx q[3];
rz(0.74362459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84393152) q[2];
sx q[2];
rz(-2.6466978) q[2];
sx q[2];
rz(2.3893791) q[2];
rz(0.080605896) q[3];
sx q[3];
rz(-1.3496496) q[3];
sx q[3];
rz(0.54752553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053631) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(-2.2930131) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(-2.2304429) q[2];
sx q[2];
rz(-1.0187889) q[2];
sx q[2];
rz(2.7762882) q[2];
rz(0.6497339) q[3];
sx q[3];
rz(-2.5841373) q[3];
sx q[3];
rz(2.0852603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
