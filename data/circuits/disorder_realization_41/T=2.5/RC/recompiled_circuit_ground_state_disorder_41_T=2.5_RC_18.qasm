OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(-1.2394152) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1151506) q[0];
sx q[0];
rz(-1.7998994) q[0];
sx q[0];
rz(-1.312448) q[0];
x q[1];
rz(0.26387604) q[2];
sx q[2];
rz(-1.7603163) q[2];
sx q[2];
rz(-0.60631982) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2762294) q[1];
sx q[1];
rz(-1.6003834) q[1];
sx q[1];
rz(-2.8113356) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69783437) q[3];
sx q[3];
rz(-0.69785944) q[3];
sx q[3];
rz(-2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.3489464) q[2];
rz(-2.3979483) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(-0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.141356) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(-2.8479688) q[0];
rz(2.1108421) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(0.34293276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4727508) q[0];
sx q[0];
rz(-2.0905417) q[0];
sx q[0];
rz(0.23730554) q[0];
rz(0.12983506) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(-0.8952199) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98951116) q[1];
sx q[1];
rz(-0.73484269) q[1];
sx q[1];
rz(2.8192725) q[1];
rz(-pi) q[2];
rz(-0.94128709) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(-0.91703892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(3.1208842) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6193806) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(-0.64747539) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-2.3562145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42703907) q[0];
sx q[0];
rz(-1.4903965) q[0];
sx q[0];
rz(-0.15958448) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9663229) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(-2.4966405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33357802) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(2.1120815) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80180577) q[3];
sx q[3];
rz(-2.2983645) q[3];
sx q[3];
rz(-1.3096732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(-2.8010098) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(-1.0964099) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88056662) q[0];
sx q[0];
rz(-1.9184904) q[0];
sx q[0];
rz(0.33190042) q[0];
rz(-pi) q[1];
rz(1.5113513) q[2];
sx q[2];
rz(-1.6013833) q[2];
sx q[2];
rz(-2.2074043) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8453522) q[1];
sx q[1];
rz(-1.3120323) q[1];
sx q[1];
rz(0.47611632) q[1];
rz(-pi) q[2];
rz(0.55620749) q[3];
sx q[3];
rz(-2.3283151) q[3];
sx q[3];
rz(-0.63268703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-0.30245936) q[0];
sx q[0];
rz(2.6746993) q[0];
rz(-0.11058552) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(2.0707524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0127481) q[0];
sx q[0];
rz(-1.332988) q[0];
sx q[0];
rz(1.1046011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51009615) q[2];
sx q[2];
rz(-2.7802711) q[2];
sx q[2];
rz(2.9053743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5683978) q[1];
sx q[1];
rz(-1.0612951) q[1];
sx q[1];
rz(2.0054818) q[1];
rz(-pi) q[2];
rz(-0.91584622) q[3];
sx q[3];
rz(-0.98036843) q[3];
sx q[3];
rz(-0.46186033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0232627) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(2.2892717) q[2];
rz(-3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0090050176) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(-1.9992794) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2347533) q[0];
sx q[0];
rz(-1.2182353) q[0];
sx q[0];
rz(-1.6749188) q[0];
rz(-pi) q[1];
rz(-3.0794607) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(-1.0203938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9276322) q[1];
sx q[1];
rz(-0.42812706) q[1];
sx q[1];
rz(-0.63251782) q[1];
rz(1.3942424) q[3];
sx q[3];
rz(-1.2640959) q[3];
sx q[3];
rz(-2.1666985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3219354) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(3.0211871) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(-2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33893809) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(1.5974207) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(2.6905751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4225578) q[0];
sx q[0];
rz(-1.8717878) q[0];
sx q[0];
rz(-1.5084356) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9092009) q[2];
sx q[2];
rz(-1.2133779) q[2];
sx q[2];
rz(0.1279624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4828795) q[1];
sx q[1];
rz(-0.45156839) q[1];
sx q[1];
rz(0.2803448) q[1];
rz(2.3312206) q[3];
sx q[3];
rz(-1.5775418) q[3];
sx q[3];
rz(2.7254259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1958127) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.7808524) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(-1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9129974) q[0];
sx q[0];
rz(-1.5902336) q[0];
sx q[0];
rz(0.0073846505) q[0];
rz(-pi) q[1];
rz(2.3191586) q[2];
sx q[2];
rz(-1.8959799) q[2];
sx q[2];
rz(-0.48516824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4256712) q[1];
sx q[1];
rz(-1.7737264) q[1];
sx q[1];
rz(0.58260609) q[1];
rz(-pi) q[2];
rz(0.13369067) q[3];
sx q[3];
rz(-2.055759) q[3];
sx q[3];
rz(-3.0669989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2989444) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(0.16743463) q[2];
rz(1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3779959) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(-1.7077712) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094291501) q[0];
sx q[0];
rz(-0.30891788) q[0];
sx q[0];
rz(-2.7500217) q[0];
rz(-pi) q[1];
rz(3.07934) q[2];
sx q[2];
rz(-0.57800284) q[2];
sx q[2];
rz(0.5731155) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3206818) q[1];
sx q[1];
rz(-1.4467941) q[1];
sx q[1];
rz(2.4895666) q[1];
rz(1.7172708) q[3];
sx q[3];
rz(-2.4333242) q[3];
sx q[3];
rz(-2.012501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47664777) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(-0.42759582) q[2];
rz(1.5852196) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(3.1184149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4172702) q[0];
sx q[0];
rz(-2.5886154) q[0];
sx q[0];
rz(-0.31674762) q[0];
x q[1];
rz(-2.542664) q[2];
sx q[2];
rz(-2.1492276) q[2];
sx q[2];
rz(0.79858649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1286298) q[1];
sx q[1];
rz(-2.9010765) q[1];
sx q[1];
rz(-0.81324767) q[1];
rz(-pi) q[2];
rz(2.3007727) q[3];
sx q[3];
rz(-1.7125907) q[3];
sx q[3];
rz(-0.26103448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(2.6507157) q[2];
rz(2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28521095) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-2.7748952) q[2];
sx q[2];
rz(-2.4187805) q[2];
sx q[2];
rz(0.4772966) q[2];
rz(1.5506977) q[3];
sx q[3];
rz(-2.1051959) q[3];
sx q[3];
rz(0.013230562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
