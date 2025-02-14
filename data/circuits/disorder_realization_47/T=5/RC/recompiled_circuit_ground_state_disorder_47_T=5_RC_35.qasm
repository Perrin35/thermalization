OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(3.0106944) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(0.17732492) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180358) q[0];
sx q[0];
rz(-2.3157488) q[0];
sx q[0];
rz(2.7808583) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2134883) q[2];
sx q[2];
rz(-0.36967247) q[2];
sx q[2];
rz(-3.0189454) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63943458) q[1];
sx q[1];
rz(-1.3356767) q[1];
sx q[1];
rz(-2.058063) q[1];
rz(-2.8683441) q[3];
sx q[3];
rz(-2.1299703) q[3];
sx q[3];
rz(2.3756053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6191972) q[2];
sx q[2];
rz(-1.088524) q[2];
sx q[2];
rz(1.8001451) q[2];
rz(-1.3522735) q[3];
sx q[3];
rz(-1.6100581) q[3];
sx q[3];
rz(1.2927829) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.751048) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(-2.6265662) q[0];
rz(2.7797049) q[1];
sx q[1];
rz(-1.7444976) q[1];
sx q[1];
rz(2.1451758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032971708) q[0];
sx q[0];
rz(-1.7713799) q[0];
sx q[0];
rz(1.8114835) q[0];
x q[1];
rz(1.1364765) q[2];
sx q[2];
rz(-1.4114209) q[2];
sx q[2];
rz(-0.61018131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39033088) q[1];
sx q[1];
rz(-1.0511025) q[1];
sx q[1];
rz(3.0355176) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89869638) q[3];
sx q[3];
rz(-0.65124159) q[3];
sx q[3];
rz(-2.8399093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7289875) q[2];
sx q[2];
rz(-2.1952486) q[2];
sx q[2];
rz(-1.6607355) q[2];
rz(-2.6442773) q[3];
sx q[3];
rz(-0.9484843) q[3];
sx q[3];
rz(2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6052674) q[0];
sx q[0];
rz(-1.241642) q[0];
sx q[0];
rz(-1.190825) q[0];
rz(0.39769998) q[1];
sx q[1];
rz(-1.560248) q[1];
sx q[1];
rz(-0.89888987) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54972285) q[0];
sx q[0];
rz(-1.6373349) q[0];
sx q[0];
rz(1.9107242) q[0];
x q[1];
rz(1.6562881) q[2];
sx q[2];
rz(-1.7608661) q[2];
sx q[2];
rz(0.23641931) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11375107) q[1];
sx q[1];
rz(-1.6154624) q[1];
sx q[1];
rz(-1.0381446) q[1];
rz(-1.8679138) q[3];
sx q[3];
rz(-2.6032902) q[3];
sx q[3];
rz(-2.3915714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1206104) q[2];
sx q[2];
rz(-1.1889428) q[2];
sx q[2];
rz(-2.9494542) q[2];
rz(2.4118679) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(-2.5016968) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925053) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(-1.8335023) q[0];
rz(-2.2970842) q[1];
sx q[1];
rz(-1.4627855) q[1];
sx q[1];
rz(0.62215296) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62646851) q[0];
sx q[0];
rz(-0.65931407) q[0];
sx q[0];
rz(0.71172442) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2842769) q[2];
sx q[2];
rz(-0.19294365) q[2];
sx q[2];
rz(-2.3920453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6250861) q[1];
sx q[1];
rz(-1.7343905) q[1];
sx q[1];
rz(-1.518599) q[1];
x q[2];
rz(-2.801311) q[3];
sx q[3];
rz(-2.8064686) q[3];
sx q[3];
rz(-1.7652896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1239803) q[2];
sx q[2];
rz(-1.3717185) q[2];
sx q[2];
rz(2.2464216) q[2];
rz(0.34919843) q[3];
sx q[3];
rz(-2.5694191) q[3];
sx q[3];
rz(2.9077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8643841) q[0];
sx q[0];
rz(-1.207749) q[0];
sx q[0];
rz(-2.5982017) q[0];
rz(0.1132938) q[1];
sx q[1];
rz(-1.3985059) q[1];
sx q[1];
rz(1.8494122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47885146) q[0];
sx q[0];
rz(-2.2355192) q[0];
sx q[0];
rz(0.49481884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0512442) q[2];
sx q[2];
rz(-0.70864973) q[2];
sx q[2];
rz(2.6330269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1317392) q[1];
sx q[1];
rz(-0.92127548) q[1];
sx q[1];
rz(-0.76503896) q[1];
x q[2];
rz(-0.73907799) q[3];
sx q[3];
rz(-1.0024286) q[3];
sx q[3];
rz(-0.06423244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8285404) q[2];
sx q[2];
rz(-1.3581759) q[2];
sx q[2];
rz(0.48946112) q[2];
rz(-0.66568565) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(2.3556975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996443) q[0];
sx q[0];
rz(-2.2388832) q[0];
sx q[0];
rz(-0.06074252) q[0];
rz(-1.2784866) q[1];
sx q[1];
rz(-0.91438952) q[1];
sx q[1];
rz(-2.1155105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64549146) q[0];
sx q[0];
rz(-2.0466986) q[0];
sx q[0];
rz(-2.2055513) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0842956) q[2];
sx q[2];
rz(-0.62556534) q[2];
sx q[2];
rz(1.4826258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85426165) q[1];
sx q[1];
rz(-1.8214487) q[1];
sx q[1];
rz(1.7118368) q[1];
rz(-2.6164651) q[3];
sx q[3];
rz(-0.59288247) q[3];
sx q[3];
rz(1.2660932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8817899) q[2];
sx q[2];
rz(-1.5851574) q[2];
sx q[2];
rz(0.078710236) q[2];
rz(1.4824661) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(0.44221529) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.215613) q[0];
sx q[0];
rz(-0.38023606) q[0];
sx q[0];
rz(-1.6967787) q[0];
rz(-2.3346057) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(-0.011946202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1769971) q[0];
sx q[0];
rz(-2.2095293) q[0];
sx q[0];
rz(-0.033319017) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4623227) q[2];
sx q[2];
rz(-1.0656992) q[2];
sx q[2];
rz(0.5031571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50004369) q[1];
sx q[1];
rz(-2.1807359) q[1];
sx q[1];
rz(-0.088492101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1720522) q[3];
sx q[3];
rz(-1.4280025) q[3];
sx q[3];
rz(-0.90313426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5936467) q[2];
sx q[2];
rz(-2.839489) q[2];
sx q[2];
rz(-2.3054403) q[2];
rz(-3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.9713822) q[0];
sx q[0];
rz(-2.2844071) q[0];
sx q[0];
rz(-0.49825391) q[0];
rz(-1.0055297) q[1];
sx q[1];
rz(-1.4509095) q[1];
sx q[1];
rz(-0.11140579) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.446602) q[0];
sx q[0];
rz(-1.1518307) q[0];
sx q[0];
rz(-0.10563157) q[0];
x q[1];
rz(-0.51014238) q[2];
sx q[2];
rz(-2.6594793) q[2];
sx q[2];
rz(3.1010266) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1131683) q[1];
sx q[1];
rz(-2.3481124) q[1];
sx q[1];
rz(-1.1641349) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7965589) q[3];
sx q[3];
rz(-1.4514231) q[3];
sx q[3];
rz(-0.81689801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19693836) q[2];
sx q[2];
rz(-1.3613746) q[2];
sx q[2];
rz(-1.7353752) q[2];
rz(-1.9841291) q[3];
sx q[3];
rz(-2.1486053) q[3];
sx q[3];
rz(1.5520613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.120753) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-1.7027759) q[0];
rz(2.2976047) q[1];
sx q[1];
rz(-1.3071209) q[1];
sx q[1];
rz(-2.9356975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015259764) q[0];
sx q[0];
rz(-1.5301989) q[0];
sx q[0];
rz(-2.5102528) q[0];
x q[1];
rz(-1.5175784) q[2];
sx q[2];
rz(-2.2578265) q[2];
sx q[2];
rz(1.1498957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83826021) q[1];
sx q[1];
rz(-0.15619677) q[1];
sx q[1];
rz(2.6172943) q[1];
x q[2];
rz(-1.1631257) q[3];
sx q[3];
rz(-2.0317997) q[3];
sx q[3];
rz(1.7940652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13731185) q[2];
sx q[2];
rz(-2.0840804) q[2];
sx q[2];
rz(1.6005969) q[2];
rz(-0.48504034) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(0.11390991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6523022) q[0];
sx q[0];
rz(-0.052209608) q[0];
sx q[0];
rz(-1.8452277) q[0];
rz(-1.0507874) q[1];
sx q[1];
rz(-1.4605582) q[1];
sx q[1];
rz(0.60660648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832336) q[0];
sx q[0];
rz(-0.16671058) q[0];
sx q[0];
rz(1.3728844) q[0];
x q[1];
rz(-1.6982618) q[2];
sx q[2];
rz(-1.0147139) q[2];
sx q[2];
rz(-0.13917222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1639599) q[1];
sx q[1];
rz(-0.31144413) q[1];
sx q[1];
rz(2.9410579) q[1];
x q[2];
rz(0.28429417) q[3];
sx q[3];
rz(-1.458567) q[3];
sx q[3];
rz(-2.1865213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7217094) q[2];
sx q[2];
rz(-1.1315283) q[2];
sx q[2];
rz(-0.69880542) q[2];
rz(-1.4612259) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(2.6644126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2496495) q[0];
sx q[0];
rz(-1.6310545) q[0];
sx q[0];
rz(-1.7153836) q[0];
rz(-0.37829788) q[1];
sx q[1];
rz(-0.62722798) q[1];
sx q[1];
rz(2.41082) q[1];
rz(-1.0964285) q[2];
sx q[2];
rz(-1.8468936) q[2];
sx q[2];
rz(1.5441128) q[2];
rz(-1.5398239) q[3];
sx q[3];
rz(-2.1270558) q[3];
sx q[3];
rz(1.1598641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
