OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(1.6311128) q[0];
rz(1.0881967) q[1];
sx q[1];
rz(-1.3047941) q[1];
sx q[1];
rz(2.3545177) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45739503) q[0];
sx q[0];
rz(-1.6826864) q[0];
sx q[0];
rz(-0.68023139) q[0];
x q[1];
rz(-2.745109) q[2];
sx q[2];
rz(-1.6210437) q[2];
sx q[2];
rz(-0.32279245) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30336994) q[1];
sx q[1];
rz(-1.359424) q[1];
sx q[1];
rz(-1.6081393) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98104279) q[3];
sx q[3];
rz(-2.1043805) q[3];
sx q[3];
rz(-0.55280815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8481019) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(0.21609406) q[2];
rz(0.50348336) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(0.066472806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1450495) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-0.37522069) q[0];
rz(0.61620617) q[1];
sx q[1];
rz(-1.0169949) q[1];
sx q[1];
rz(-2.9527051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045563) q[0];
sx q[0];
rz(-2.3052373) q[0];
sx q[0];
rz(-1.0009074) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1470931) q[2];
sx q[2];
rz(-0.81953555) q[2];
sx q[2];
rz(0.34162765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0068854) q[1];
sx q[1];
rz(-1.7733679) q[1];
sx q[1];
rz(0.88500513) q[1];
rz(-0.62633212) q[3];
sx q[3];
rz(-2.5181558) q[3];
sx q[3];
rz(1.745519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45794332) q[2];
sx q[2];
rz(-1.318149) q[2];
sx q[2];
rz(1.9981492) q[2];
rz(-0.6360561) q[3];
sx q[3];
rz(-3.0885185) q[3];
sx q[3];
rz(-1.599357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77883333) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(-1.1676316) q[0];
rz(0.12655839) q[1];
sx q[1];
rz(-1.626222) q[1];
sx q[1];
rz(-0.57356858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9605947) q[0];
sx q[0];
rz(-2.3947236) q[0];
sx q[0];
rz(1.7043714) q[0];
rz(-2.7748108) q[2];
sx q[2];
rz(-1.8213118) q[2];
sx q[2];
rz(2.9652852) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8178278) q[1];
sx q[1];
rz(-2.7080302) q[1];
sx q[1];
rz(1.4563515) q[1];
rz(-pi) q[2];
rz(-2.6010547) q[3];
sx q[3];
rz(-2.4063568) q[3];
sx q[3];
rz(1.0395887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48587376) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(-1.5029079) q[2];
rz(1.8466628) q[3];
sx q[3];
rz(-2.86627) q[3];
sx q[3];
rz(0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5824222) q[0];
sx q[0];
rz(-2.9952413) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(-2.9811663) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(-0.697335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9458185) q[0];
sx q[0];
rz(-1.8911288) q[0];
sx q[0];
rz(0.26024241) q[0];
rz(-1.2658046) q[2];
sx q[2];
rz(-2.2100497) q[2];
sx q[2];
rz(1.5196619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26665822) q[1];
sx q[1];
rz(-1.8968079) q[1];
sx q[1];
rz(-0.74649109) q[1];
rz(-pi) q[2];
rz(-0.13015161) q[3];
sx q[3];
rz(-1.5846545) q[3];
sx q[3];
rz(-2.7133301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9852898) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-1.1227597) q[2];
rz(1.5924234) q[3];
sx q[3];
rz(-1.4607818) q[3];
sx q[3];
rz(2.7482225) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3093981) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(-1.1291809) q[0];
rz(-2.26217) q[1];
sx q[1];
rz(-1.3112473) q[1];
sx q[1];
rz(2.5130491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.578772) q[0];
sx q[0];
rz(-1.6878753) q[0];
sx q[0];
rz(1.6952008) q[0];
x q[1];
rz(2.1206003) q[2];
sx q[2];
rz(-2.1537184) q[2];
sx q[2];
rz(-1.1958808) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0776808) q[1];
sx q[1];
rz(-1.6110629) q[1];
sx q[1];
rz(2.1230554) q[1];
x q[2];
rz(1.699963) q[3];
sx q[3];
rz(-2.3979275) q[3];
sx q[3];
rz(2.5088333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-0.76801378) q[2];
sx q[2];
rz(-1.3735695) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-1.0765358) q[3];
sx q[3];
rz(3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3304928) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(-2.9789441) q[0];
rz(-1.2341518) q[1];
sx q[1];
rz(-1.994543) q[1];
sx q[1];
rz(-0.92076921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4402078) q[0];
sx q[0];
rz(-0.21344859) q[0];
sx q[0];
rz(-1.0769072) q[0];
rz(2.2170794) q[2];
sx q[2];
rz(-2.0343896) q[2];
sx q[2];
rz(1.6789503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9879976) q[1];
sx q[1];
rz(-0.82057602) q[1];
sx q[1];
rz(-1.4819891) q[1];
rz(-pi) q[2];
rz(0.1699888) q[3];
sx q[3];
rz(-0.78944474) q[3];
sx q[3];
rz(1.1657749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84772253) q[2];
sx q[2];
rz(-0.72124481) q[2];
sx q[2];
rz(-3.0976307) q[2];
rz(1.4219159) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(-3.1066084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.1304355) q[0];
sx q[0];
rz(-2.3373117) q[0];
sx q[0];
rz(-3.0479808) q[0];
rz(-2.1162107) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(-2.3536918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073613361) q[0];
sx q[0];
rz(-2.1755621) q[0];
sx q[0];
rz(-0.66901274) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56886832) q[2];
sx q[2];
rz(-1.7143357) q[2];
sx q[2];
rz(1.6678866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.390139) q[1];
sx q[1];
rz(-2.1691704) q[1];
sx q[1];
rz(1.5214439) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1206081) q[3];
sx q[3];
rz(-1.4292307) q[3];
sx q[3];
rz(1.0673351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0660144) q[2];
sx q[2];
rz(-0.54328537) q[2];
sx q[2];
rz(-2.1755966) q[2];
rz(2.994359) q[3];
sx q[3];
rz(-1.6817663) q[3];
sx q[3];
rz(0.60432965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(1.2526441) q[0];
rz(0.057627536) q[1];
sx q[1];
rz(-1.7689972) q[1];
sx q[1];
rz(-0.39841121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1712043) q[0];
sx q[0];
rz(-2.1792767) q[0];
sx q[0];
rz(1.7826016) q[0];
rz(-pi) q[1];
rz(-0.49968991) q[2];
sx q[2];
rz(-1.4213287) q[2];
sx q[2];
rz(-2.2563344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7207253) q[1];
sx q[1];
rz(-1.5667672) q[1];
sx q[1];
rz(1.3708877) q[1];
x q[2];
rz(-0.54315217) q[3];
sx q[3];
rz(-1.6622322) q[3];
sx q[3];
rz(-2.9989797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5742699) q[2];
sx q[2];
rz(-1.2719354) q[2];
sx q[2];
rz(0.92180139) q[2];
rz(-2.4343893) q[3];
sx q[3];
rz(-1.0450398) q[3];
sx q[3];
rz(-0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1929753) q[0];
sx q[0];
rz(-0.12117584) q[0];
sx q[0];
rz(-2.2344672) q[0];
rz(2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9655749) q[0];
sx q[0];
rz(-1.2621573) q[0];
sx q[0];
rz(1.465331) q[0];
rz(-pi) q[1];
rz(-2.3263518) q[2];
sx q[2];
rz(-2.4324907) q[2];
sx q[2];
rz(0.39295024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.231016) q[1];
sx q[1];
rz(-1.8635995) q[1];
sx q[1];
rz(-0.80033719) q[1];
rz(-2.4828047) q[3];
sx q[3];
rz(-1.4843936) q[3];
sx q[3];
rz(0.98443809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85764641) q[2];
sx q[2];
rz(-1.593677) q[2];
sx q[2];
rz(2.471981) q[2];
rz(-2.5745463) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(-1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6794716) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(-0.8959499) q[0];
rz(-3.1013007) q[1];
sx q[1];
rz(-0.94172421) q[1];
sx q[1];
rz(1.5774073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98151796) q[0];
sx q[0];
rz(-1.5943564) q[0];
sx q[0];
rz(0.0047163211) q[0];
x q[1];
rz(2.0032004) q[2];
sx q[2];
rz(-2.2096425) q[2];
sx q[2];
rz(-0.8919208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9312694) q[1];
sx q[1];
rz(-2.0788553) q[1];
sx q[1];
rz(-0.47048752) q[1];
x q[2];
rz(2.778042) q[3];
sx q[3];
rz(-2.145316) q[3];
sx q[3];
rz(-2.1839104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3080052) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(-1.3099111) q[2];
rz(-1.1854019) q[3];
sx q[3];
rz(-0.88806051) q[3];
sx q[3];
rz(1.5230702) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014521) q[0];
sx q[0];
rz(-1.3012713) q[0];
sx q[0];
rz(-0.95300994) q[0];
rz(0.078527191) q[1];
sx q[1];
rz(-1.086906) q[1];
sx q[1];
rz(-0.79798098) q[1];
rz(-0.25987249) q[2];
sx q[2];
rz(-0.91283471) q[2];
sx q[2];
rz(1.2531493) q[2];
rz(2.4918741) q[3];
sx q[3];
rz(-0.25593723) q[3];
sx q[3];
rz(-1.2091936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
