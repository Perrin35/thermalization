OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2406727) q[0];
sx q[0];
rz(-2.9958041) q[0];
sx q[0];
rz(-1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9300728) q[0];
sx q[0];
rz(-2.2317108) q[0];
sx q[0];
rz(-1.5050864) q[0];
rz(0.47122987) q[2];
sx q[2];
rz(-1.8719851) q[2];
sx q[2];
rz(0.1567086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92857526) q[1];
sx q[1];
rz(-1.8452106) q[1];
sx q[1];
rz(-1.3540512) q[1];
rz(-pi) q[2];
rz(0.59009976) q[3];
sx q[3];
rz(-1.3076772) q[3];
sx q[3];
rz(2.2424169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(0.8068223) q[2];
rz(0.72871366) q[3];
sx q[3];
rz(-0.7586793) q[3];
sx q[3];
rz(1.9369283) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(-0.30701315) q[0];
rz(-1.2767876) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(2.0397287) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80230882) q[0];
sx q[0];
rz(-1.3972939) q[0];
sx q[0];
rz(2.1306778) q[0];
rz(-2.8028433) q[2];
sx q[2];
rz(-2.951606) q[2];
sx q[2];
rz(-0.39443406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15436048) q[1];
sx q[1];
rz(-1.7061966) q[1];
sx q[1];
rz(-1.3920648) q[1];
x q[2];
rz(2.3086433) q[3];
sx q[3];
rz(-2.2112339) q[3];
sx q[3];
rz(0.73329496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1133984) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(3.0451575) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.6343445) q[3];
sx q[3];
rz(-0.69795394) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(-1.6290889) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(2.8834744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494875) q[0];
sx q[0];
rz(-1.4319812) q[0];
sx q[0];
rz(-1.3347244) q[0];
x q[1];
rz(-1.6506972) q[2];
sx q[2];
rz(-1.6703484) q[2];
sx q[2];
rz(-0.71074206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8363179) q[1];
sx q[1];
rz(-1.5757676) q[1];
sx q[1];
rz(3.1286376) q[1];
rz(1.3985004) q[3];
sx q[3];
rz(-1.3329643) q[3];
sx q[3];
rz(-0.83816499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.021598024) q[2];
sx q[2];
rz(-1.1979411) q[2];
sx q[2];
rz(0.032133948) q[2];
rz(-2.8739127) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699919) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(3.0573523) q[0];
rz(0.035942297) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(2.8004004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223136) q[0];
sx q[0];
rz(-1.3094762) q[0];
sx q[0];
rz(2.4909291) q[0];
x q[1];
rz(0.3772328) q[2];
sx q[2];
rz(-2.5209941) q[2];
sx q[2];
rz(-2.008174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0964716) q[1];
sx q[1];
rz(-2.0092138) q[1];
sx q[1];
rz(-2.3271266) q[1];
rz(-0.59935948) q[3];
sx q[3];
rz(-1.1903569) q[3];
sx q[3];
rz(2.3886556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(-1.535447) q[2];
rz(0.26220194) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(0.95482993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3839805) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(2.0641932) q[0];
rz(-0.39240882) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(2.1108625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2868256) q[0];
sx q[0];
rz(-1.0374026) q[0];
sx q[0];
rz(1.0817762) q[0];
rz(-pi) q[1];
rz(-1.0723128) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(-1.0637525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64068078) q[1];
sx q[1];
rz(-1.7927153) q[1];
sx q[1];
rz(2.9725916) q[1];
rz(-pi) q[2];
rz(-0.34319539) q[3];
sx q[3];
rz(-1.0846018) q[3];
sx q[3];
rz(-0.23517683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82421676) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(-2.2281036) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(-0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(-1.6656026) q[0];
rz(-0.39235517) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(2.5501693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77350835) q[0];
sx q[0];
rz(-2.4249833) q[0];
sx q[0];
rz(2.1518097) q[0];
x q[1];
rz(1.3900969) q[2];
sx q[2];
rz(-1.8417674) q[2];
sx q[2];
rz(-2.9621552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6063664) q[1];
sx q[1];
rz(-0.51035817) q[1];
sx q[1];
rz(-1.7061069) q[1];
rz(-0.37844946) q[3];
sx q[3];
rz(-0.90988084) q[3];
sx q[3];
rz(-1.701041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-0.66606194) q[2];
sx q[2];
rz(0.57233125) q[2];
rz(-0.22219292) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(2.4967994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-0.40827665) q[0];
rz(2.4093742) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(-0.29762038) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2985122) q[0];
sx q[0];
rz(-1.1491927) q[0];
sx q[0];
rz(1.1888191) q[0];
rz(-pi) q[1];
rz(-1.6082798) q[2];
sx q[2];
rz(-0.49502326) q[2];
sx q[2];
rz(1.9062454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0875594) q[1];
sx q[1];
rz(-0.39759025) q[1];
sx q[1];
rz(2.1000186) q[1];
rz(-1.7209531) q[3];
sx q[3];
rz(-0.87267002) q[3];
sx q[3];
rz(2.0719178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(2.4453898) q[2];
rz(-2.2512186) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625229) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-0.21275511) q[0];
rz(-0.46956024) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(0.75417095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.629166) q[0];
sx q[0];
rz(-0.34722024) q[0];
sx q[0];
rz(-0.99033333) q[0];
x q[1];
rz(0.32525678) q[2];
sx q[2];
rz(-1.3768679) q[2];
sx q[2];
rz(0.12530279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27867815) q[1];
sx q[1];
rz(-2.2154741) q[1];
sx q[1];
rz(-1.4281322) q[1];
x q[2];
rz(-2.4741094) q[3];
sx q[3];
rz(-2.1060447) q[3];
sx q[3];
rz(0.014674295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.135005) q[2];
sx q[2];
rz(-0.90234119) q[2];
sx q[2];
rz(0.78224409) q[2];
rz(-1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-0.35000354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049147216) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(0.96442047) q[0];
rz(-1.8655221) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(1.5020348) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0014711) q[0];
sx q[0];
rz(-1.3706511) q[0];
sx q[0];
rz(-1.3339304) q[0];
rz(-pi) q[1];
rz(-1.1310857) q[2];
sx q[2];
rz(-2.2520116) q[2];
sx q[2];
rz(2.0063673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.561475) q[1];
sx q[1];
rz(-2.0881542) q[1];
sx q[1];
rz(-3.1151732) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.067599452) q[3];
sx q[3];
rz(-1.4238333) q[3];
sx q[3];
rz(0.33892469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8883349) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(0.9453195) q[2];
rz(2.3156598) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(2.6065684) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0655521) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(2.3384576) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.6604661) q[1];
sx q[1];
rz(0.28958431) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591551) q[0];
sx q[0];
rz(-1.6763601) q[0];
sx q[0];
rz(-3.1292874) q[0];
x q[1];
rz(1.2425735) q[2];
sx q[2];
rz(-1.2699091) q[2];
sx q[2];
rz(0.93133486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9152569) q[1];
sx q[1];
rz(-2.7136972) q[1];
sx q[1];
rz(2.1386599) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8554162) q[3];
sx q[3];
rz(-1.2499785) q[3];
sx q[3];
rz(-0.30485728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-2.181459) q[2];
rz(-2.5586186) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(2.2768903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.690602) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(0.040738978) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(0.30754752) q[2];
sx q[2];
rz(-1.1196305) q[2];
sx q[2];
rz(1.011871) q[2];
rz(1.4996281) q[3];
sx q[3];
rz(-1.6215848) q[3];
sx q[3];
rz(3.0231089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
