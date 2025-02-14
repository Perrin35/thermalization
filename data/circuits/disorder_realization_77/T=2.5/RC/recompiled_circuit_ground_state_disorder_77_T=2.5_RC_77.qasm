OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9935432) q[0];
sx q[0];
rz(-2.7033959) q[0];
sx q[0];
rz(-0.73102695) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(-2.5875523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95931388) q[0];
sx q[0];
rz(-2.4327299) q[0];
sx q[0];
rz(-0.76226632) q[0];
rz(2.0494387) q[2];
sx q[2];
rz(-1.3502099) q[2];
sx q[2];
rz(-1.5696862) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3202127) q[1];
sx q[1];
rz(-2.2578635) q[1];
sx q[1];
rz(-0.043619556) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1393806) q[3];
sx q[3];
rz(-2.3750288) q[3];
sx q[3];
rz(-2.1994906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76250166) q[2];
sx q[2];
rz(-1.5662301) q[2];
sx q[2];
rz(-2.9762034) q[2];
rz(-0.23073828) q[3];
sx q[3];
rz(-2.8985891) q[3];
sx q[3];
rz(2.2802584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6623401) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(-1.3586556) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-0.75444573) q[1];
sx q[1];
rz(0.84017909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0884893) q[0];
sx q[0];
rz(-1.8376956) q[0];
sx q[0];
rz(-2.0574112) q[0];
rz(-0.86960881) q[2];
sx q[2];
rz(-1.879862) q[2];
sx q[2];
rz(2.3562507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11586861) q[1];
sx q[1];
rz(-2.3626973) q[1];
sx q[1];
rz(0.42515691) q[1];
rz(-pi) q[2];
rz(1.4235557) q[3];
sx q[3];
rz(-0.81063945) q[3];
sx q[3];
rz(-0.37732201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75300616) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(0.93442717) q[2];
rz(-1.2855444) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(0.88750315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11338209) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(-1.5895948) q[0];
rz(1.7734843) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(2.0210463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5277953) q[0];
sx q[0];
rz(-2.5585735) q[0];
sx q[0];
rz(-2.1979419) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8875966) q[2];
sx q[2];
rz(-1.1448556) q[2];
sx q[2];
rz(-2.6563702) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9160794) q[1];
sx q[1];
rz(-0.54359964) q[1];
sx q[1];
rz(-2.9617642) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1638223) q[3];
sx q[3];
rz(-1.6864136) q[3];
sx q[3];
rz(1.8979929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0418732) q[2];
sx q[2];
rz(-0.21328829) q[2];
sx q[2];
rz(2.8495157) q[2];
rz(-1.9000351) q[3];
sx q[3];
rz(-1.7009267) q[3];
sx q[3];
rz(2.6888473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9860155) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(2.7098932) q[0];
rz(-2.9875634) q[1];
sx q[1];
rz(-1.5715716) q[1];
sx q[1];
rz(-1.0607176) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78846473) q[0];
sx q[0];
rz(-1.9110838) q[0];
sx q[0];
rz(-2.8022604) q[0];
rz(-pi) q[1];
rz(-0.17825408) q[2];
sx q[2];
rz(-1.3298103) q[2];
sx q[2];
rz(0.18496938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1338443) q[1];
sx q[1];
rz(-0.40002003) q[1];
sx q[1];
rz(2.176214) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4134365) q[3];
sx q[3];
rz(-1.3829136) q[3];
sx q[3];
rz(1.431501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74956191) q[2];
sx q[2];
rz(-1.531484) q[2];
sx q[2];
rz(-2.5677666) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.6618238) q[3];
sx q[3];
rz(-0.81006947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.5350128) q[0];
sx q[0];
rz(-0.25595328) q[0];
sx q[0];
rz(2.8888597) q[0];
rz(-0.17929721) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(1.4089233) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104467) q[0];
sx q[0];
rz(-2.8928693) q[0];
sx q[0];
rz(-1.876271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9581355) q[2];
sx q[2];
rz(-1.9692067) q[2];
sx q[2];
rz(0.18926316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45917967) q[1];
sx q[1];
rz(-0.79303592) q[1];
sx q[1];
rz(0.17211087) q[1];
x q[2];
rz(0.5770251) q[3];
sx q[3];
rz(-0.60363704) q[3];
sx q[3];
rz(1.9848422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0130284) q[2];
sx q[2];
rz(-1.7594254) q[2];
sx q[2];
rz(1.9405091) q[2];
rz(-0.80738336) q[3];
sx q[3];
rz(-1.3599334) q[3];
sx q[3];
rz(2.0090296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(3.0627237) q[0];
sx q[0];
rz(-1.4987334) q[0];
sx q[0];
rz(-2.3003182) q[0];
rz(-1.9367283) q[1];
sx q[1];
rz(-1.8236022) q[1];
sx q[1];
rz(-1.0296317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.556419) q[0];
sx q[0];
rz(-1.5812567) q[0];
sx q[0];
rz(2.0812698) q[0];
rz(0.51635833) q[2];
sx q[2];
rz(-1.8022924) q[2];
sx q[2];
rz(1.4524492) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8358546) q[1];
sx q[1];
rz(-1.5630013) q[1];
sx q[1];
rz(-2.9628997) q[1];
x q[2];
rz(-2.5978686) q[3];
sx q[3];
rz(-0.35055509) q[3];
sx q[3];
rz(-1.2389099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.010217696) q[2];
sx q[2];
rz(-0.75675941) q[2];
sx q[2];
rz(-2.6436464) q[2];
rz(-0.52305269) q[3];
sx q[3];
rz(-1.7224576) q[3];
sx q[3];
rz(0.12652346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.842857) q[0];
sx q[0];
rz(-0.13164483) q[0];
sx q[0];
rz(0.81197062) q[0];
rz(-0.89093351) q[1];
sx q[1];
rz(-1.1633326) q[1];
sx q[1];
rz(0.99937159) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97702128) q[0];
sx q[0];
rz(-1.1081105) q[0];
sx q[0];
rz(-2.0624731) q[0];
rz(-pi) q[1];
rz(-0.096316263) q[2];
sx q[2];
rz(-1.673398) q[2];
sx q[2];
rz(1.4136537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4564086) q[1];
sx q[1];
rz(-1.0899836) q[1];
sx q[1];
rz(-0.64688869) q[1];
x q[2];
rz(1.9346775) q[3];
sx q[3];
rz(-0.69857222) q[3];
sx q[3];
rz(-1.1766124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.823395) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(1.8434175) q[2];
rz(3.1031109) q[3];
sx q[3];
rz(-2.0352071) q[3];
sx q[3];
rz(1.5255671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0865974) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(0.31563345) q[0];
rz(-1.9075958) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(-1.2589781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4988853) q[0];
sx q[0];
rz(-2.2872637) q[0];
sx q[0];
rz(-1.1906719) q[0];
rz(0.76541111) q[2];
sx q[2];
rz(-1.4452626) q[2];
sx q[2];
rz(-1.5620934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96871829) q[1];
sx q[1];
rz(-0.58191381) q[1];
sx q[1];
rz(0.7919148) q[1];
rz(-pi) q[2];
x q[2];
rz(0.057825967) q[3];
sx q[3];
rz(-0.10198051) q[3];
sx q[3];
rz(-1.1241871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2235609) q[2];
sx q[2];
rz(-2.3980902) q[2];
sx q[2];
rz(-0.027916748) q[2];
rz(-2.3935086) q[3];
sx q[3];
rz(-1.4192162) q[3];
sx q[3];
rz(-2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41392031) q[0];
sx q[0];
rz(-0.99069178) q[0];
sx q[0];
rz(-0.45898166) q[0];
rz(-2.0052295) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(0.21839011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39762625) q[0];
sx q[0];
rz(-1.5874001) q[0];
sx q[0];
rz(-3.1182609) q[0];
rz(3.0074429) q[2];
sx q[2];
rz(-1.9512259) q[2];
sx q[2];
rz(2.5594437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9105007) q[1];
sx q[1];
rz(-2.0606961) q[1];
sx q[1];
rz(1.7375768) q[1];
rz(-pi) q[2];
rz(-1.9588542) q[3];
sx q[3];
rz(-1.4493296) q[3];
sx q[3];
rz(1.1282008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2961262) q[2];
sx q[2];
rz(-2.8870236) q[2];
sx q[2];
rz(2.0938342) q[2];
rz(0.74538499) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(1.6489702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202241) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-0.34287232) q[0];
rz(2.951237) q[1];
sx q[1];
rz(-2.1326667) q[1];
sx q[1];
rz(-0.51317936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032470908) q[0];
sx q[0];
rz(-2.4801804) q[0];
sx q[0];
rz(1.3107857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97810271) q[2];
sx q[2];
rz(-0.58916559) q[2];
sx q[2];
rz(1.9798673) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.733963) q[1];
sx q[1];
rz(-1.9701951) q[1];
sx q[1];
rz(-2.6435889) q[1];
x q[2];
rz(-1.8796092) q[3];
sx q[3];
rz(-0.43918375) q[3];
sx q[3];
rz(-2.9787946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7036983) q[2];
sx q[2];
rz(-0.63592211) q[2];
sx q[2];
rz(0.61100125) q[2];
rz(-0.25887394) q[3];
sx q[3];
rz(-1.2114108) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3872869) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(-0.98987956) q[1];
sx q[1];
rz(-1.9342593) q[1];
sx q[1];
rz(0.63607279) q[1];
rz(-1.5237332) q[2];
sx q[2];
rz(-0.18885352) q[2];
sx q[2];
rz(0.39299008) q[2];
rz(2.4425735) q[3];
sx q[3];
rz(-1.5200281) q[3];
sx q[3];
rz(2.716223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
