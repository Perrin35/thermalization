OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(-1.3949431) q[0];
sx q[0];
rz(-3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.476386) q[0];
sx q[0];
rz(-1.6205661) q[0];
sx q[0];
rz(2.9988891) q[0];
x q[1];
rz(-2.7156419) q[2];
sx q[2];
rz(-0.89086878) q[2];
sx q[2];
rz(-1.8217063) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.98382681) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(1.8450792) q[1];
rz(-pi) q[2];
rz(-0.28489057) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(0.26998664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.3886064) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(-1.0999854) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.3756479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(0.21582614) q[0];
x q[1];
rz(-1.3865115) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(-2.266303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.980629) q[1];
sx q[1];
rz(-1.7501608) q[1];
sx q[1];
rz(0.37100002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0166753) q[3];
sx q[3];
rz(-1.8801873) q[3];
sx q[3];
rz(0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21330825) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(2.0080163) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.2352357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2429758) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(2.5234733) q[0];
rz(-pi) q[1];
rz(0.58358242) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(-0.32354087) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5879844) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(-1.8334465) q[1];
rz(2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(-0.94961387) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(-3.0217357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(3.0981307) q[0];
rz(-0.36074952) q[2];
sx q[2];
rz(-1.4348239) q[2];
sx q[2];
rz(-0.13775682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90384342) q[1];
sx q[1];
rz(-1.3386054) q[1];
sx q[1];
rz(3.1234427) q[1];
x q[2];
rz(0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2729028) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(-0.23916434) q[0];
rz(0.24057062) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(2.0672928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6606635) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(1.6515345) q[1];
rz(2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.6606768) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-3.1408177) q[2];
sx q[2];
rz(-1.1256071) q[2];
sx q[2];
rz(-2.0489401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.1793544) q[1];
x q[2];
rz(0.89208608) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(-0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(-2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(-2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.3247103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0990471) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(3.0142473) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84682805) q[2];
sx q[2];
rz(-0.70038751) q[2];
sx q[2];
rz(-2.2221153) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.054338) q[1];
sx q[1];
rz(-1.9458276) q[1];
sx q[1];
rz(2.7021728) q[1];
x q[2];
rz(0.38447325) q[3];
sx q[3];
rz(-0.95015804) q[3];
sx q[3];
rz(-1.4409161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(-2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136918) q[0];
sx q[0];
rz(-2.0223589) q[0];
sx q[0];
rz(-2.5149462) q[0];
rz(-0.82215674) q[2];
sx q[2];
rz(-2.2347054) q[2];
sx q[2];
rz(-2.8234931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76028484) q[1];
sx q[1];
rz(-0.21796255) q[1];
sx q[1];
rz(2.1078307) q[1];
x q[2];
rz(1.7361705) q[3];
sx q[3];
rz(-0.91455063) q[3];
sx q[3];
rz(2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(0.28265488) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219014) q[0];
sx q[0];
rz(-1.1150517) q[0];
sx q[0];
rz(-0.68463188) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33090584) q[2];
sx q[2];
rz(-1.7754284) q[2];
sx q[2];
rz(0.76588878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54284755) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(0.41136841) q[1];
x q[2];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(-3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0081351) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(-1.0159147) q[0];
x q[1];
rz(-2.290906) q[2];
sx q[2];
rz(-0.88973532) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1055923) q[1];
sx q[1];
rz(-2.3899374) q[1];
sx q[1];
rz(2.4113301) q[1];
rz(-pi) q[2];
rz(0.82716771) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(-0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(-1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.66424673) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(1.3689465) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];