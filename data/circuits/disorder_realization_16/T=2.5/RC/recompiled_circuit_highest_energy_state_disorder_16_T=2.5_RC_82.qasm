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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(2.733736) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(-1.8899625) q[1];
sx q[1];
rz(-1.6012021) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8046869) q[0];
sx q[0];
rz(-1.1338455) q[0];
sx q[0];
rz(3.0880465) q[0];
rz(-pi) q[1];
rz(-2.8694334) q[2];
sx q[2];
rz(-2.6725884) q[2];
sx q[2];
rz(-0.16527612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1931128) q[1];
sx q[1];
rz(-2.6958637) q[1];
sx q[1];
rz(-0.50909252) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82117053) q[3];
sx q[3];
rz(-1.0945012) q[3];
sx q[3];
rz(1.1007635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5539598) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(2.9239192) q[2];
rz(2.830128) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(1.1658839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199453) q[0];
sx q[0];
rz(-0.2758652) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(2.3852589) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(-1.5914894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645351) q[0];
sx q[0];
rz(-1.4486396) q[0];
sx q[0];
rz(2.4273901) q[0];
rz(1.1728376) q[2];
sx q[2];
rz(-1.6670456) q[2];
sx q[2];
rz(-2.725696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8988425) q[1];
sx q[1];
rz(-1.6926584) q[1];
sx q[1];
rz(2.980848) q[1];
rz(-pi) q[2];
rz(0.85568537) q[3];
sx q[3];
rz(-2.4484903) q[3];
sx q[3];
rz(2.2732796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99593607) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(-0.23925979) q[2];
rz(-1.650882) q[3];
sx q[3];
rz(-1.8473293) q[3];
sx q[3];
rz(1.2283121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2759129) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-3.1269585) q[0];
rz(-2.2485661) q[1];
sx q[1];
rz(-0.97517401) q[1];
sx q[1];
rz(0.21569529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0365643) q[0];
sx q[0];
rz(-1.6220548) q[0];
sx q[0];
rz(-3.1369484) q[0];
rz(-pi) q[1];
rz(-1.9142308) q[2];
sx q[2];
rz(-0.94007713) q[2];
sx q[2];
rz(-1.7523927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5657422) q[1];
sx q[1];
rz(-1.8025959) q[1];
sx q[1];
rz(1.7017468) q[1];
rz(-pi) q[2];
rz(-2.2641597) q[3];
sx q[3];
rz(-2.139353) q[3];
sx q[3];
rz(-1.7753778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0821685) q[2];
sx q[2];
rz(-1.9614204) q[2];
sx q[2];
rz(-2.1480985) q[2];
rz(-0.70837402) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(-3.0181001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3941536) q[0];
sx q[0];
rz(-2.6469632) q[0];
sx q[0];
rz(0.69394773) q[0];
rz(-0.28678647) q[1];
sx q[1];
rz(-1.5319805) q[1];
sx q[1];
rz(1.087711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74831731) q[0];
sx q[0];
rz(-1.6946185) q[0];
sx q[0];
rz(1.0626777) q[0];
rz(1.0469646) q[2];
sx q[2];
rz(-1.0817384) q[2];
sx q[2];
rz(-0.2917052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36632167) q[1];
sx q[1];
rz(-1.2676116) q[1];
sx q[1];
rz(-2.1949185) q[1];
x q[2];
rz(1.5252211) q[3];
sx q[3];
rz(-0.85698313) q[3];
sx q[3];
rz(2.8136106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9168758) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(-1.7737596) q[2];
rz(1.9035089) q[3];
sx q[3];
rz(-1.5725458) q[3];
sx q[3];
rz(2.9253173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5991768) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(-2.5237778) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-2.8384659) q[1];
sx q[1];
rz(-0.88315001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1533494) q[0];
sx q[0];
rz(-0.81479615) q[0];
sx q[0];
rz(1.5606461) q[0];
rz(-pi) q[1];
rz(2.8034535) q[2];
sx q[2];
rz(-0.45881841) q[2];
sx q[2];
rz(0.79566075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8700734) q[1];
sx q[1];
rz(-0.49527676) q[1];
sx q[1];
rz(1.4325159) q[1];
rz(2.7297121) q[3];
sx q[3];
rz(-1.05384) q[3];
sx q[3];
rz(0.98983562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21657476) q[2];
sx q[2];
rz(-0.25187945) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(-1.2733634) q[3];
sx q[3];
rz(-1.4342156) q[3];
sx q[3];
rz(-0.98289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044416044) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(-0.62028766) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(1.1995859) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9387377) q[0];
sx q[0];
rz(-1.5498383) q[0];
sx q[0];
rz(-1.6771132) q[0];
rz(-2.1753005) q[2];
sx q[2];
rz(-0.25560616) q[2];
sx q[2];
rz(-0.57744712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5127677) q[1];
sx q[1];
rz(-0.56124748) q[1];
sx q[1];
rz(0.36833771) q[1];
rz(-pi) q[2];
rz(0.32034822) q[3];
sx q[3];
rz(-0.36133063) q[3];
sx q[3];
rz(-0.02905127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8951796) q[2];
sx q[2];
rz(-1.8975039) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(-2.8504168) q[3];
sx q[3];
rz(-2.5305735) q[3];
sx q[3];
rz(0.62180716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596844) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(2.9631462) q[0];
rz(-2.2715691) q[1];
sx q[1];
rz(-2.2585637) q[1];
sx q[1];
rz(-0.80284405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69235301) q[0];
sx q[0];
rz(-0.62004161) q[0];
sx q[0];
rz(-0.043278261) q[0];
rz(0.58391352) q[2];
sx q[2];
rz(-1.3656934) q[2];
sx q[2];
rz(-2.5896304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72910129) q[1];
sx q[1];
rz(-1.8961398) q[1];
sx q[1];
rz(-0.71112432) q[1];
x q[2];
rz(-1.9196627) q[3];
sx q[3];
rz(-2.4838944) q[3];
sx q[3];
rz(-1.5546834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7650083) q[2];
sx q[2];
rz(-1.8439801) q[2];
sx q[2];
rz(-2.2479748) q[2];
rz(-1.5997959) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1641418) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(-0.2151016) q[0];
rz(-2.2970301) q[1];
sx q[1];
rz(-1.8212049) q[1];
sx q[1];
rz(-0.79427687) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85958033) q[0];
sx q[0];
rz(-1.5826125) q[0];
sx q[0];
rz(0.026547308) q[0];
rz(-pi) q[1];
rz(1.2634981) q[2];
sx q[2];
rz(-2.0086882) q[2];
sx q[2];
rz(-1.4222933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6805598) q[1];
sx q[1];
rz(-0.54745142) q[1];
sx q[1];
rz(1.7458818) q[1];
rz(1.2450386) q[3];
sx q[3];
rz(-0.81607258) q[3];
sx q[3];
rz(-1.2118424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.27791417) q[2];
sx q[2];
rz(-1.2028376) q[2];
sx q[2];
rz(-1.4097144) q[2];
rz(-0.75285161) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(-2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.79733217) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(-2.1687188) q[0];
rz(-1.5219888) q[1];
sx q[1];
rz(-1.7771959) q[1];
sx q[1];
rz(0.63046986) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75062932) q[0];
sx q[0];
rz(-1.4328151) q[0];
sx q[0];
rz(-0.18355455) q[0];
rz(1.5233598) q[2];
sx q[2];
rz(-2.9844797) q[2];
sx q[2];
rz(-0.6873129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8119295) q[1];
sx q[1];
rz(-2.3721937) q[1];
sx q[1];
rz(-3.0698464) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9679873) q[3];
sx q[3];
rz(-1.6328873) q[3];
sx q[3];
rz(0.28854077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.015532739) q[2];
sx q[2];
rz(-1.1522747) q[2];
sx q[2];
rz(-1.8401592) q[2];
rz(-1.556501) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(2.7263156) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179317) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(2.6483722) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.247765) q[1];
sx q[1];
rz(-2.7769322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64283687) q[0];
sx q[0];
rz(-2.7793573) q[0];
sx q[0];
rz(1.9693768) q[0];
rz(1.0027867) q[2];
sx q[2];
rz(-0.73633654) q[2];
sx q[2];
rz(-1.9970837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2233878) q[1];
sx q[1];
rz(-2.3406696) q[1];
sx q[1];
rz(0.078322874) q[1];
x q[2];
rz(-1.6935891) q[3];
sx q[3];
rz(-2.8229575) q[3];
sx q[3];
rz(2.2991869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6898592) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(-0.9922007) q[2];
rz(-0.36710468) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(-0.67830694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.622396) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(0.92338152) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(-0.95465701) q[2];
sx q[2];
rz(-2.8980394) q[2];
sx q[2];
rz(0.63799636) q[2];
rz(-2.4319875) q[3];
sx q[3];
rz(-2.3313739) q[3];
sx q[3];
rz(-0.45562638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
