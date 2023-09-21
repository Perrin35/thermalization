OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5231397) q[0];
sx q[0];
rz(-1.8123527) q[0];
sx q[0];
rz(-3.0735344) q[0];
rz(-pi) q[1];
rz(-2.152359) q[2];
sx q[2];
rz(-0.49058149) q[2];
sx q[2];
rz(1.4197592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0393385) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(2.0233872) q[1];
x q[2];
rz(-1.8688698) q[3];
sx q[3];
rz(-2.094659) q[3];
sx q[3];
rz(0.03447547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0704437) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-2.0106993) q[0];
x q[1];
rz(-1.8663835) q[2];
sx q[2];
rz(-0.81381961) q[2];
sx q[2];
rz(2.5676167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21287316) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(-2.5961155) q[1];
rz(-pi) q[2];
rz(-1.133425) q[3];
sx q[3];
rz(-1.8873894) q[3];
sx q[3];
rz(-0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(2.6142696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2505328) q[0];
sx q[0];
rz(-0.30936229) q[0];
sx q[0];
rz(-3.1191349) q[0];
rz(-pi) q[1];
rz(1.2450561) q[2];
sx q[2];
rz(-0.92391652) q[2];
sx q[2];
rz(-1.5145472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(0.36348344) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35964386) q[3];
sx q[3];
rz(-0.76562866) q[3];
sx q[3];
rz(-1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-2.3153268) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6986615) q[0];
sx q[0];
rz(-0.94416617) q[0];
sx q[0];
rz(2.4311964) q[0];
rz(-1.516953) q[2];
sx q[2];
rz(-0.87140897) q[2];
sx q[2];
rz(2.2211423) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5845881) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(-1.9083379) q[1];
x q[2];
rz(2.5097333) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782398) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(2.1468303) q[0];
rz(-pi) q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(-2.9885459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51296556) q[1];
sx q[1];
rz(-1.8497397) q[1];
sx q[1];
rz(0.97475027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17554749) q[3];
sx q[3];
rz(-1.5834337) q[3];
sx q[3];
rz(2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(0.0897952) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3072309) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(2.2063072) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1744556) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(-2.9784163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.60702) q[1];
sx q[1];
rz(-1.2828476) q[1];
sx q[1];
rz(-1.7400017) q[1];
x q[2];
rz(-2.0635707) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.1759261) q[2];
rz(3.0684493) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(2.6426962) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38567625) q[0];
sx q[0];
rz(-1.809281) q[0];
sx q[0];
rz(0.083906108) q[0];
rz(-1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(1.6183311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(0.1715626) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8764898) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508142) q[0];
sx q[0];
rz(-2.6065126) q[0];
sx q[0];
rz(-2.4909766) q[0];
rz(-2.5152399) q[2];
sx q[2];
rz(-0.92712958) q[2];
sx q[2];
rz(0.49382892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9951524) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(-0.050683024) q[1];
rz(0.96945854) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.2649149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3919968) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(-1.0438265) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0090946) q[2];
sx q[2];
rz(-1.0602789) q[2];
sx q[2];
rz(0.27062182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1688655) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(1.6070073) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7183652) q[3];
sx q[3];
rz(-1.3864577) q[3];
sx q[3];
rz(-1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.4642749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20477644) q[0];
sx q[0];
rz(-1.4602772) q[0];
sx q[0];
rz(-0.21893455) q[0];
x q[1];
rz(1.9015058) q[2];
sx q[2];
rz(-0.61243528) q[2];
sx q[2];
rz(-1.2182747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6357395) q[1];
sx q[1];
rz(-0.74843279) q[1];
sx q[1];
rz(-1.1670508) q[1];
rz(-pi) q[2];
rz(-3.0562923) q[3];
sx q[3];
rz(-1.7461667) q[3];
sx q[3];
rz(-2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(2.5901569) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(2.5682156) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
