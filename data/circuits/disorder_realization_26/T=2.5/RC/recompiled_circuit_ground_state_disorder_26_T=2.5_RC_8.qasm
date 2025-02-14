OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5698009) q[0];
sx q[0];
rz(4.0507841) q[0];
sx q[0];
rz(7.9180766) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(6.7758898) q[1];
sx q[1];
rz(13.125782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052033612) q[0];
sx q[0];
rz(-1.8843643) q[0];
sx q[0];
rz(-1.7455717) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.962553) q[2];
sx q[2];
rz(-1.9843352) q[2];
sx q[2];
rz(-0.79038436) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8570429) q[1];
sx q[1];
rz(-2.8618097) q[1];
sx q[1];
rz(-1.5955052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1640638) q[3];
sx q[3];
rz(-1.1479605) q[3];
sx q[3];
rz(-0.55227597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59014615) q[2];
sx q[2];
rz(-0.85942736) q[2];
sx q[2];
rz(2.0476511) q[2];
rz(-1.9048196) q[3];
sx q[3];
rz(-1.9779132) q[3];
sx q[3];
rz(2.8804603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7673489) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(-0.71794024) q[0];
rz(-1.9473437) q[1];
sx q[1];
rz(-2.7337044) q[1];
sx q[1];
rz(-1.6494707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46280086) q[0];
sx q[0];
rz(-2.959842) q[0];
sx q[0];
rz(-0.75580718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70121558) q[2];
sx q[2];
rz(-0.55169902) q[2];
sx q[2];
rz(3.015446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1528984) q[1];
sx q[1];
rz(-2.0190658) q[1];
sx q[1];
rz(2.6397735) q[1];
rz(0.2160875) q[3];
sx q[3];
rz(-2.3041445) q[3];
sx q[3];
rz(2.2186389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0481999) q[2];
sx q[2];
rz(-0.77646774) q[2];
sx q[2];
rz(0.36273599) q[2];
rz(2.2705966) q[3];
sx q[3];
rz(-1.3716776) q[3];
sx q[3];
rz(1.8179651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094630346) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(-0.56075019) q[0];
rz(-0.97460711) q[1];
sx q[1];
rz(-2.2798996) q[1];
sx q[1];
rz(2.9240756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95642521) q[0];
sx q[0];
rz(-2.5568058) q[0];
sx q[0];
rz(-0.034791481) q[0];
x q[1];
rz(-1.0571805) q[2];
sx q[2];
rz(-1.660778) q[2];
sx q[2];
rz(-1.0124026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8629093) q[1];
sx q[1];
rz(-0.28537286) q[1];
sx q[1];
rz(-1.1527658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6337745) q[3];
sx q[3];
rz(-0.96270409) q[3];
sx q[3];
rz(-2.5127895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64323032) q[2];
sx q[2];
rz(-0.78780323) q[2];
sx q[2];
rz(-0.95334774) q[2];
rz(-2.0638454) q[3];
sx q[3];
rz(-1.6809623) q[3];
sx q[3];
rz(2.6673711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0039907) q[0];
sx q[0];
rz(-2.9676262) q[0];
sx q[0];
rz(-0.91651383) q[0];
rz(-0.42463955) q[1];
sx q[1];
rz(-1.759513) q[1];
sx q[1];
rz(2.0133846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5461499) q[0];
sx q[0];
rz(-0.91888035) q[0];
sx q[0];
rz(2.2007887) q[0];
rz(-pi) q[1];
rz(-1.7776971) q[2];
sx q[2];
rz(-0.82447663) q[2];
sx q[2];
rz(-1.4109703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2353246) q[1];
sx q[1];
rz(-2.5282994) q[1];
sx q[1];
rz(2.6032531) q[1];
rz(-pi) q[2];
rz(0.22596328) q[3];
sx q[3];
rz(-0.98528242) q[3];
sx q[3];
rz(-0.33285347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8757223) q[2];
sx q[2];
rz(-1.46572) q[2];
sx q[2];
rz(1.1992559) q[2];
rz(-1.9843598) q[3];
sx q[3];
rz(-2.3374989) q[3];
sx q[3];
rz(-2.6257302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1542094) q[0];
sx q[0];
rz(-0.043954285) q[0];
sx q[0];
rz(0.92535812) q[0];
rz(-0.56272733) q[1];
sx q[1];
rz(-1.0147164) q[1];
sx q[1];
rz(-2.7664807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35532668) q[0];
sx q[0];
rz(-2.5857175) q[0];
sx q[0];
rz(-2.9340266) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0882224) q[2];
sx q[2];
rz(-2.0056021) q[2];
sx q[2];
rz(-2.82719) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8567849) q[1];
sx q[1];
rz(-1.7397659) q[1];
sx q[1];
rz(-1.595975) q[1];
rz(-1.1697606) q[3];
sx q[3];
rz(-1.9936863) q[3];
sx q[3];
rz(2.324375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23842266) q[2];
sx q[2];
rz(-2.2245202) q[2];
sx q[2];
rz(0.49752107) q[2];
rz(0.43429747) q[3];
sx q[3];
rz(-1.4562166) q[3];
sx q[3];
rz(2.1527877) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76181805) q[0];
sx q[0];
rz(-2.2803545) q[0];
sx q[0];
rz(0.27989835) q[0];
rz(0.99413904) q[1];
sx q[1];
rz(-0.86052624) q[1];
sx q[1];
rz(-0.57847374) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8168479) q[0];
sx q[0];
rz(-1.9473796) q[0];
sx q[0];
rz(2.5927932) q[0];
rz(0.45409338) q[2];
sx q[2];
rz(-0.75579772) q[2];
sx q[2];
rz(0.86802717) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38505618) q[1];
sx q[1];
rz(-2.2900683) q[1];
sx q[1];
rz(-1.3100708) q[1];
x q[2];
rz(-1.469645) q[3];
sx q[3];
rz(-1.8378075) q[3];
sx q[3];
rz(-1.8570516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0824739) q[2];
sx q[2];
rz(-0.59933496) q[2];
sx q[2];
rz(-2.2590051) q[2];
rz(-2.5146218) q[3];
sx q[3];
rz(-2.4866703) q[3];
sx q[3];
rz(-0.89491189) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6249228) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(0.00061568419) q[0];
rz(0.90283886) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(-0.045086233) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1182528) q[0];
sx q[0];
rz(-2.059142) q[0];
sx q[0];
rz(1.4834542) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3995724) q[2];
sx q[2];
rz(-1.7270403) q[2];
sx q[2];
rz(2.7458423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2704986) q[1];
sx q[1];
rz(-1.2648858) q[1];
sx q[1];
rz(2.3853777) q[1];
rz(-pi) q[2];
rz(-0.99211971) q[3];
sx q[3];
rz(-2.1014155) q[3];
sx q[3];
rz(-1.19095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8467329) q[2];
sx q[2];
rz(-0.555987) q[2];
sx q[2];
rz(-0.9168469) q[2];
rz(-2.2461241) q[3];
sx q[3];
rz(-1.6296891) q[3];
sx q[3];
rz(1.2112613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13670066) q[0];
sx q[0];
rz(-3.0558375) q[0];
sx q[0];
rz(2.6766747) q[0];
rz(-1.9526941) q[1];
sx q[1];
rz(-1.7949972) q[1];
sx q[1];
rz(2.7704923) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3329778) q[0];
sx q[0];
rz(-0.59691256) q[0];
sx q[0];
rz(2.3521168) q[0];
rz(0.82000374) q[2];
sx q[2];
rz(-1.5981673) q[2];
sx q[2];
rz(-3.0153831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5587204) q[1];
sx q[1];
rz(-1.1338738) q[1];
sx q[1];
rz(-0.82804273) q[1];
rz(-0.63432548) q[3];
sx q[3];
rz(-0.46503572) q[3];
sx q[3];
rz(0.1022235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66172877) q[2];
sx q[2];
rz(-0.78833818) q[2];
sx q[2];
rz(2.7313477) q[2];
rz(1.6346301) q[3];
sx q[3];
rz(-2.7362636) q[3];
sx q[3];
rz(0.043225616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0931382) q[0];
sx q[0];
rz(-0.52670902) q[0];
sx q[0];
rz(2.9658537) q[0];
rz(2.3161092) q[1];
sx q[1];
rz(-1.261542) q[1];
sx q[1];
rz(-2.3186191) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12619647) q[0];
sx q[0];
rz(-1.0312197) q[0];
sx q[0];
rz(-0.71643512) q[0];
rz(-2.1954567) q[2];
sx q[2];
rz(-2.6310134) q[2];
sx q[2];
rz(-0.389314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8205679) q[1];
sx q[1];
rz(-1.2234283) q[1];
sx q[1];
rz(0.72183164) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9814616) q[3];
sx q[3];
rz(-1.6045827) q[3];
sx q[3];
rz(0.50212348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8601816) q[2];
sx q[2];
rz(-1.9965636) q[2];
sx q[2];
rz(0.54110503) q[2];
rz(2.7152854) q[3];
sx q[3];
rz(-0.46190327) q[3];
sx q[3];
rz(2.2947252) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4538039) q[0];
sx q[0];
rz(-0.80950469) q[0];
sx q[0];
rz(-0.33568207) q[0];
rz(0.022739284) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(-2.4679599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27972066) q[0];
sx q[0];
rz(-2.4422944) q[0];
sx q[0];
rz(0.55681105) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.672253) q[2];
sx q[2];
rz(-1.1038269) q[2];
sx q[2];
rz(1.4604258) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1991829) q[1];
sx q[1];
rz(-1.7734151) q[1];
sx q[1];
rz(-0.97572787) q[1];
x q[2];
rz(1.7106179) q[3];
sx q[3];
rz(-1.2396024) q[3];
sx q[3];
rz(-3.0013623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.444904) q[2];
sx q[2];
rz(-0.99385571) q[2];
sx q[2];
rz(1.0059086) q[2];
rz(2.3575947) q[3];
sx q[3];
rz(-2.8204212) q[3];
sx q[3];
rz(-1.706749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1226596) q[0];
sx q[0];
rz(-1.4228595) q[0];
sx q[0];
rz(2.0056437) q[0];
rz(-0.38416531) q[1];
sx q[1];
rz(-1.9728248) q[1];
sx q[1];
rz(-1.493175) q[1];
rz(0.76569414) q[2];
sx q[2];
rz(-2.8509344) q[2];
sx q[2];
rz(0.40727587) q[2];
rz(2.0060282) q[3];
sx q[3];
rz(-0.76261997) q[3];
sx q[3];
rz(-0.95773209) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
