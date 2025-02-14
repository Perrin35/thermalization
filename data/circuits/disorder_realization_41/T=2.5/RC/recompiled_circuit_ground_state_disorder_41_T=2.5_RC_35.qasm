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
rz(2.4089101) q[1];
sx q[1];
rz(10.664193) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(-0.23668134) q[0];
x q[1];
rz(0.26387604) q[2];
sx q[2];
rz(-1.7603163) q[2];
sx q[2];
rz(-0.60631982) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2762294) q[1];
sx q[1];
rz(-1.5412093) q[1];
sx q[1];
rz(0.33025708) q[1];
x q[2];
rz(2.4437583) q[3];
sx q[3];
rz(-0.69785944) q[3];
sx q[3];
rz(2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4702845) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(1.7926463) q[2];
rz(-0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(0.17717895) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002366) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(2.8479688) q[0];
rz(1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(2.7986599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2156649) q[0];
sx q[0];
rz(-0.56680381) q[0];
sx q[0];
rz(-1.180992) q[0];
x q[1];
rz(1.8937102) q[2];
sx q[2];
rz(-0.39034778) q[2];
sx q[2];
rz(-0.54803145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.575042) q[1];
sx q[1];
rz(-2.2600265) q[1];
sx q[1];
rz(-1.8495967) q[1];
rz(-2.8906451) q[3];
sx q[3];
rz(-0.95635575) q[3];
sx q[3];
rz(0.80048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0127516) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-0.7592321) q[2];
rz(-0.020708474) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(0.52221209) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(3.1371327) q[0];
rz(0.64747539) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42703907) q[0];
sx q[0];
rz(-1.4903965) q[0];
sx q[0];
rz(-0.15958448) q[0];
x q[1];
rz(2.9663229) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(-0.64495211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2783616) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(0.079915483) q[1];
rz(2.3397869) q[3];
sx q[3];
rz(-2.2983645) q[3];
sx q[3];
rz(-1.8319195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(-1.7837589) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(-2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99712813) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(1.0964099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344656) q[0];
sx q[0];
rz(-1.2594481) q[0];
sx q[0];
rz(1.9368571) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6302413) q[2];
sx q[2];
rz(-1.6013833) q[2];
sx q[2];
rz(2.2074043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29624048) q[1];
sx q[1];
rz(-1.3120323) q[1];
sx q[1];
rz(2.6654763) q[1];
rz(-pi) q[2];
rz(-1.0616333) q[3];
sx q[3];
rz(-0.90583767) q[3];
sx q[3];
rz(0.10275118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61802822) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(0.63557449) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(-1.0708403) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3239646) q[0];
sx q[0];
rz(-1.1187129) q[0];
sx q[0];
rz(-0.26497193) q[0];
rz(1.3883287) q[2];
sx q[2];
rz(-1.8844205) q[2];
sx q[2];
rz(-2.8387866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5731949) q[1];
sx q[1];
rz(-1.0612951) q[1];
sx q[1];
rz(1.1361109) q[1];
x q[2];
rz(-2.2257464) q[3];
sx q[3];
rz(-2.1612242) q[3];
sx q[3];
rz(2.6797323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(2.2892717) q[2];
rz(-3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(1.8527385) q[0];
rz(1.6291078) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.9992794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5286251) q[0];
sx q[0];
rz(-2.7745947) q[0];
sx q[0];
rz(-0.27530833) q[0];
rz(-pi) q[1];
rz(0.062131957) q[2];
sx q[2];
rz(-1.3634063) q[2];
sx q[2];
rz(2.1211989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2494275) q[1];
sx q[1];
rz(-1.9122423) q[1];
sx q[1];
rz(1.83431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3942424) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(-2.1666985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.8737996) q[2];
sx q[2];
rz(-0.12040559) q[2];
rz(1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.45101756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2149725) q[0];
sx q[0];
rz(-0.30719137) q[0];
sx q[0];
rz(-2.9434669) q[0];
rz(-pi) q[1];
rz(2.9092009) q[2];
sx q[2];
rz(-1.9282148) q[2];
sx q[2];
rz(0.1279624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.34141065) q[1];
sx q[1];
rz(-1.4497611) q[1];
sx q[1];
rz(2.7054663) q[1];
rz(-pi) q[2];
rz(0.0093098442) q[3];
sx q[3];
rz(-0.81039372) q[3];
sx q[3];
rz(-1.9805465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1958127) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(-1.3607402) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(2.4355167) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(2.2854038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5917011) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(2.3191586) q[2];
sx q[2];
rz(-1.8959799) q[2];
sx q[2];
rz(2.6564244) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7159214) q[1];
sx q[1];
rz(-1.3678663) q[1];
sx q[1];
rz(0.58260609) q[1];
x q[2];
rz(-2.059465) q[3];
sx q[3];
rz(-1.4525982) q[3];
sx q[3];
rz(1.4335872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(2.974158) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.7077712) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3146798) q[0];
sx q[0];
rz(-1.2859435) q[0];
sx q[0];
rz(1.4495984) q[0];
rz(-0.062252684) q[2];
sx q[2];
rz(-0.57800284) q[2];
sx q[2];
rz(0.5731155) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65574291) q[1];
sx q[1];
rz(-2.2169737) q[1];
sx q[1];
rz(-1.7263401) q[1];
rz(3.0172273) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(1.82064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47664777) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(0.023177711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231004) q[0];
sx q[0];
rz(-1.4064624) q[0];
sx q[0];
rz(0.53043764) q[0];
x q[1];
rz(0.59892861) q[2];
sx q[2];
rz(-0.99236503) q[2];
sx q[2];
rz(-0.79858649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35649037) q[1];
sx q[1];
rz(-1.39686) q[1];
sx q[1];
rz(-2.9746303) q[1];
rz(-pi) q[2];
rz(2.9523207) q[3];
sx q[3];
rz(-2.2918275) q[3];
sx q[3];
rz(1.4355575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8699441) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(-0.49087697) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8563817) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-2.8648227) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-0.68885531) q[2];
sx q[2];
rz(-1.3313455) q[2];
sx q[2];
rz(-0.81305885) q[2];
rz(-1.5908949) q[3];
sx q[3];
rz(-2.1051959) q[3];
sx q[3];
rz(0.013230562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
