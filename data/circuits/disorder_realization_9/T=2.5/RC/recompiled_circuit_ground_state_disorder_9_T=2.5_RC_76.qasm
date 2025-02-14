OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0754492) q[0];
sx q[0];
rz(5.2392449) q[0];
sx q[0];
rz(6.2934937) q[0];
rz(-0.97996867) q[1];
sx q[1];
rz(4.5865321) q[1];
sx q[1];
rz(10.001339) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.071237) q[0];
sx q[0];
rz(-0.33215392) q[0];
sx q[0];
rz(-2.7624056) q[0];
rz(-0.61277436) q[2];
sx q[2];
rz(-0.98721993) q[2];
sx q[2];
rz(0.60223168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.986728) q[1];
sx q[1];
rz(-2.0948334) q[1];
sx q[1];
rz(-1.3681332) q[1];
rz(-pi) q[2];
rz(1.3862085) q[3];
sx q[3];
rz(-2.3999676) q[3];
sx q[3];
rz(0.0024992873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6171241) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(-1.2855533) q[2];
rz(1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(1.0514528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.224757) q[0];
sx q[0];
rz(0.24573627) q[0];
rz(2.083678) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(2.8071075) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6545821) q[0];
sx q[0];
rz(-1.0791057) q[0];
sx q[0];
rz(-1.0500748) q[0];
rz(2.7203015) q[2];
sx q[2];
rz(-0.80688804) q[2];
sx q[2];
rz(-0.56299201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2776371) q[1];
sx q[1];
rz(-2.5647128) q[1];
sx q[1];
rz(-2.3052892) q[1];
rz(1.2859551) q[3];
sx q[3];
rz(-2.1495616) q[3];
sx q[3];
rz(-2.9841686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1841396) q[2];
sx q[2];
rz(-2.2440971) q[2];
sx q[2];
rz(1.3857566) q[2];
rz(2.1615248) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053213) q[0];
sx q[0];
rz(-2.1846117) q[0];
sx q[0];
rz(1.3901688) q[0];
rz(-2.7888489) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(0.62044755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.432029) q[0];
sx q[0];
rz(-1.480446) q[0];
sx q[0];
rz(1.4960775) q[0];
rz(-pi) q[1];
rz(0.10398002) q[2];
sx q[2];
rz(-0.98538387) q[2];
sx q[2];
rz(-0.81987655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57261833) q[1];
sx q[1];
rz(-0.44562045) q[1];
sx q[1];
rz(2.3072412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6246715) q[3];
sx q[3];
rz(-2.1955588) q[3];
sx q[3];
rz(-2.3464835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1227526) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(-1.8943465) q[2];
rz(-0.16417575) q[3];
sx q[3];
rz(-1.7069867) q[3];
sx q[3];
rz(-2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71488798) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(-1.2870652) q[0];
rz(-2.5413068) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716887) q[0];
sx q[0];
rz(-1.6780417) q[0];
sx q[0];
rz(1.5041758) q[0];
rz(-1.6389334) q[2];
sx q[2];
rz(-1.1916416) q[2];
sx q[2];
rz(-0.013163002) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0541573) q[1];
sx q[1];
rz(-1.3174926) q[1];
sx q[1];
rz(1.3852055) q[1];
rz(-pi) q[2];
rz(1.2490396) q[3];
sx q[3];
rz(-1.1776973) q[3];
sx q[3];
rz(2.2571795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6220182) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(-2.3681417) q[2];
rz(-1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(-2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89161038) q[0];
sx q[0];
rz(-1.0050499) q[0];
sx q[0];
rz(2.6494359) q[0];
rz(2.6026169) q[1];
sx q[1];
rz(-1.69918) q[1];
sx q[1];
rz(-1.0999365) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4770319) q[0];
sx q[0];
rz(-2.7072621) q[0];
sx q[0];
rz(-0.6566027) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5861584) q[2];
sx q[2];
rz(-2.5298497) q[2];
sx q[2];
rz(1.7811012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7485986) q[1];
sx q[1];
rz(-1.6637319) q[1];
sx q[1];
rz(-2.3307021) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6585369) q[3];
sx q[3];
rz(-1.277458) q[3];
sx q[3];
rz(1.4774069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8010572) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(-0.56279969) q[2];
rz(-0.69989145) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7971147) q[0];
sx q[0];
rz(-2.4176702) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(-1.2190602) q[1];
sx q[1];
rz(-1.9025758) q[1];
sx q[1];
rz(0.452279) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.723123) q[0];
sx q[0];
rz(-2.4112066) q[0];
sx q[0];
rz(2.1900221) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3269749) q[2];
sx q[2];
rz(-0.53485188) q[2];
sx q[2];
rz(1.2592821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86231316) q[1];
sx q[1];
rz(-1.24839) q[1];
sx q[1];
rz(-3.1161948) q[1];
rz(1.4283435) q[3];
sx q[3];
rz(-2.1948084) q[3];
sx q[3];
rz(2.595866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(0.13709489) q[2];
rz(1.5591722) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(-0.77663511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95098507) q[0];
sx q[0];
rz(-2.3793716) q[0];
sx q[0];
rz(1.5964339) q[0];
rz(2.6243788) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(2.1545765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1008489) q[0];
sx q[0];
rz(-0.13253875) q[0];
sx q[0];
rz(-0.1610456) q[0];
x q[1];
rz(2.0275063) q[2];
sx q[2];
rz(-2.4825077) q[2];
sx q[2];
rz(0.50625077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1068808) q[1];
sx q[1];
rz(-1.8000523) q[1];
sx q[1];
rz(-0.27193141) q[1];
x q[2];
rz(-1.3201041) q[3];
sx q[3];
rz(-0.32816988) q[3];
sx q[3];
rz(-1.8712957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.78642693) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(-3.1332341) q[2];
rz(2.9110294) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(1.6647388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01345988) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(1.531456) q[0];
rz(1.6186591) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(-1.5787554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0268036) q[0];
sx q[0];
rz(-2.1675637) q[0];
sx q[0];
rz(-1.7472302) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7152856) q[2];
sx q[2];
rz(-1.943271) q[2];
sx q[2];
rz(-2.5691751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4304428) q[1];
sx q[1];
rz(-2.4949269) q[1];
sx q[1];
rz(-0.73902392) q[1];
rz(2.4736797) q[3];
sx q[3];
rz(-1.3156943) q[3];
sx q[3];
rz(0.15984838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3030887) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(-3.1040891) q[2];
rz(0.23669067) q[3];
sx q[3];
rz(-2.762837) q[3];
sx q[3];
rz(-1.2934575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26838747) q[0];
sx q[0];
rz(-0.78998843) q[0];
sx q[0];
rz(1.0733806) q[0];
rz(-0.34183303) q[1];
sx q[1];
rz(-2.5624202) q[1];
sx q[1];
rz(-0.62172186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1960012) q[0];
sx q[0];
rz(-2.669793) q[0];
sx q[0];
rz(-2.5425172) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4474117) q[2];
sx q[2];
rz(-2.2272791) q[2];
sx q[2];
rz(-0.8379762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1248884) q[1];
sx q[1];
rz(-1.0811662) q[1];
sx q[1];
rz(-1.917065) q[1];
x q[2];
rz(0.14071669) q[3];
sx q[3];
rz(-1.6575812) q[3];
sx q[3];
rz(2.0279036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5198034) q[2];
sx q[2];
rz(-2.9534464) q[2];
sx q[2];
rz(0.56358799) q[2];
rz(1.311519) q[3];
sx q[3];
rz(-1.9026285) q[3];
sx q[3];
rz(-2.3254976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1821063) q[0];
sx q[0];
rz(-2.2725548) q[0];
sx q[0];
rz(-0.8748138) q[0];
rz(1.733571) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(-2.5061238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4548774) q[0];
sx q[0];
rz(-1.2736763) q[0];
sx q[0];
rz(-0.77154205) q[0];
x q[1];
rz(2.8671625) q[2];
sx q[2];
rz(-1.9727025) q[2];
sx q[2];
rz(0.91584554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2848031) q[1];
sx q[1];
rz(-1.9724047) q[1];
sx q[1];
rz(2.1339586) q[1];
rz(-2.3395679) q[3];
sx q[3];
rz(-0.72154183) q[3];
sx q[3];
rz(2.3156543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94329876) q[2];
sx q[2];
rz(-2.0512927) q[2];
sx q[2];
rz(2.3401006) q[2];
rz(-1.21579) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(-0.041778684) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53032482) q[0];
sx q[0];
rz(-1.7014736) q[0];
sx q[0];
rz(0.3076719) q[0];
rz(1.2904185) q[1];
sx q[1];
rz(-0.94480521) q[1];
sx q[1];
rz(0.31029846) q[1];
rz(1.7062832) q[2];
sx q[2];
rz(-1.3904962) q[2];
sx q[2];
rz(2.4615859) q[2];
rz(2.484816) q[3];
sx q[3];
rz(-1.6776424) q[3];
sx q[3];
rz(1.3959891) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
