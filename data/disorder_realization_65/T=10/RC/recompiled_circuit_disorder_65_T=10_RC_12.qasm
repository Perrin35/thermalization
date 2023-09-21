OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91905347) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(-0.28199621) q[0];
rz(-2.9825748) q[2];
sx q[2];
rz(-2.3454925) q[2];
sx q[2];
rz(3.0167992) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97779146) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(1.6519283) q[1];
x q[2];
rz(-2.0243458) q[3];
sx q[3];
rz(-2.0587066) q[3];
sx q[3];
rz(-0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.5011903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879724) q[0];
sx q[0];
rz(-3.0164218) q[0];
sx q[0];
rz(0.20010389) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1019727) q[2];
sx q[2];
rz(-1.8509682) q[2];
sx q[2];
rz(-0.75129347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.048763927) q[1];
sx q[1];
rz(-1.3991014) q[1];
sx q[1];
rz(0.6789536) q[1];
x q[2];
rz(-1.9618481) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(-0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(0.93018053) q[0];
x q[1];
rz(-2.2321646) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(-0.94240377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55206628) q[1];
sx q[1];
rz(-2.8376736) q[1];
sx q[1];
rz(-2.9986831) q[1];
x q[2];
rz(3.010473) q[3];
sx q[3];
rz(-0.91851202) q[3];
sx q[3];
rz(0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.9925041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968397) q[0];
sx q[0];
rz(-1.4595932) q[0];
sx q[0];
rz(-0.24223321) q[0];
x q[1];
rz(0.43196584) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(0.84749046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035605343) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(-1.993202) q[1];
x q[2];
rz(-2.3991688) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(-2.1267668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(-1.8956005) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(-1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(0.88422424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0852172) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(-0.18116829) q[0];
x q[1];
rz(-1.1341368) q[2];
sx q[2];
rz(-1.1982802) q[2];
sx q[2];
rz(0.30581747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.65518889) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(-0.43227355) q[1];
rz(-pi) q[2];
rz(2.103881) q[3];
sx q[3];
rz(-1.8928877) q[3];
sx q[3];
rz(1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-2.172519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3450619) q[0];
sx q[0];
rz(-0.11632761) q[0];
sx q[0];
rz(1.9274812) q[0];
x q[1];
rz(-0.078209608) q[2];
sx q[2];
rz(-2.4259896) q[2];
sx q[2];
rz(-2.2323334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25989306) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(2.0395181) q[1];
rz(-pi) q[2];
x q[2];
rz(8*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(-1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51450729) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(-0.43760854) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.4356027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317658) q[0];
sx q[0];
rz(-1.0477715) q[0];
sx q[0];
rz(-2.564389) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(2.6798623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0457397) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(2.2687001) q[1];
x q[2];
rz(-1.0649879) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(1.3169469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.7763604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2030905) q[0];
sx q[0];
rz(-1.7068958) q[0];
sx q[0];
rz(-3.0654728) q[0];
rz(-1.6985967) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(2.4335361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(-2.609842) q[1];
rz(2.7216464) q[3];
sx q[3];
rz(-1.7405602) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(2.8839135) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(-0.92393595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(-2.5328013) q[0];
rz(1.5588435) q[2];
sx q[2];
rz(-2.0593615) q[2];
sx q[2];
rz(-1.2221297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20272217) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(-2.2321781) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4531035) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(-1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.7473934) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(-0.0069847981) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(-1.5241745) q[0];
x q[1];
rz(3.0496809) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(3.0222169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2529777) q[1];
sx q[1];
rz(-1.1974317) q[1];
sx q[1];
rz(0.044815973) q[1];
rz(-2.2323654) q[3];
sx q[3];
rz(-1.070676) q[3];
sx q[3];
rz(-1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.9772298) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-3.0315728) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(-2.8244143) q[3];
sx q[3];
rz(-2.2472897) q[3];
sx q[3];
rz(-0.030285611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
