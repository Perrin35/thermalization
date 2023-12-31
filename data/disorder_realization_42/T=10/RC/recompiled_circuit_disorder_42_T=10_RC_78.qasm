OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(-2.7789814) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(-0.4904823) q[1];
sx q[1];
rz(-0.34163707) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1385961) q[0];
sx q[0];
rz(-1.9435104) q[0];
sx q[0];
rz(-1.4390104) q[0];
rz(-pi) q[1];
rz(2.5506637) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(-0.17609827) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8056148) q[1];
sx q[1];
rz(-1.059638) q[1];
sx q[1];
rz(2.7989945) q[1];
rz(-pi) q[2];
rz(0.2239286) q[3];
sx q[3];
rz(-1.9609946) q[3];
sx q[3];
rz(1.3003365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5477156) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(-0.95648742) q[2];
rz(-2.9603413) q[3];
sx q[3];
rz(-0.58877188) q[3];
sx q[3];
rz(1.6199002) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0797121) q[0];
sx q[0];
rz(-3.0759838) q[0];
sx q[0];
rz(-1.1454426) q[0];
rz(-1.068813) q[1];
sx q[1];
rz(-2.309598) q[1];
sx q[1];
rz(0.72584814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1918068) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(1.9907065) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63273301) q[2];
sx q[2];
rz(-0.69721141) q[2];
sx q[2];
rz(2.7998507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5686381) q[1];
sx q[1];
rz(-1.8796088) q[1];
sx q[1];
rz(-3.1152578) q[1];
rz(-pi) q[2];
x q[2];
rz(1.782136) q[3];
sx q[3];
rz(-0.95892116) q[3];
sx q[3];
rz(0.26821995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-0.99728161) q[2];
rz(-1.4128134) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-0.93878448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7230351) q[0];
sx q[0];
rz(-2.1756873) q[0];
sx q[0];
rz(1.1285271) q[0];
rz(-1.3035125) q[1];
sx q[1];
rz(-0.729527) q[1];
sx q[1];
rz(-2.2361141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008464) q[0];
sx q[0];
rz(-1.8119537) q[0];
sx q[0];
rz(1.3748037) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1378023) q[2];
sx q[2];
rz(-0.86647063) q[2];
sx q[2];
rz(-1.1906884) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.803373) q[1];
sx q[1];
rz(-1.6446911) q[1];
sx q[1];
rz(2.282663) q[1];
rz(-pi) q[2];
rz(0.043025322) q[3];
sx q[3];
rz(-0.85775162) q[3];
sx q[3];
rz(-2.713664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7604312) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(2.1543489) q[2];
rz(-0.045771249) q[3];
sx q[3];
rz(-1.2922623) q[3];
sx q[3];
rz(0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10953294) q[0];
sx q[0];
rz(-2.4592472) q[0];
sx q[0];
rz(-3.0155244) q[0];
rz(-0.18445045) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2950738) q[0];
sx q[0];
rz(-0.82058883) q[0];
sx q[0];
rz(-1.2942765) q[0];
x q[1];
rz(1.3912958) q[2];
sx q[2];
rz(-2.0677535) q[2];
sx q[2];
rz(-0.77979445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0427525) q[1];
sx q[1];
rz(-0.71430695) q[1];
sx q[1];
rz(2.7318098) q[1];
x q[2];
rz(-0.64658029) q[3];
sx q[3];
rz(-2.1255891) q[3];
sx q[3];
rz(3.0559029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8644774) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-0.63151044) q[2];
rz(2.946092) q[3];
sx q[3];
rz(-1.2771527) q[3];
sx q[3];
rz(2.7235532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773848) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-0.46708333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5669117) q[0];
sx q[0];
rz(-1.3710638) q[0];
sx q[0];
rz(-1.1916222) q[0];
rz(-pi) q[1];
rz(2.5555243) q[2];
sx q[2];
rz(-1.7826826) q[2];
sx q[2];
rz(-0.52473247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5284164) q[1];
sx q[1];
rz(-1.7883374) q[1];
sx q[1];
rz(2.0030641) q[1];
rz(-pi) q[2];
rz(1.3777556) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(0.65929268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(2.6494027) q[2];
rz(2.8594033) q[3];
sx q[3];
rz(-1.7084833) q[3];
sx q[3];
rz(-1.5658437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2714587) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(1.4280691) q[1];
sx q[1];
rz(-2.1336887) q[1];
sx q[1];
rz(1.6451947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7487885) q[0];
sx q[0];
rz(-1.3747871) q[0];
sx q[0];
rz(-2.1013573) q[0];
rz(-0.24106579) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(-2.966553) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.306957) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(-1.6389695) q[1];
rz(0.71293998) q[3];
sx q[3];
rz(-1.5804277) q[3];
sx q[3];
rz(-2.9603017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66203403) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(1.6408287) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-0.57290572) q[0];
sx q[0];
rz(-1.0923882) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(-0.011627442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250925) q[0];
sx q[0];
rz(-1.9511127) q[0];
sx q[0];
rz(-1.8561383) q[0];
rz(-pi) q[1];
rz(-3.1161357) q[2];
sx q[2];
rz(-1.5989132) q[2];
sx q[2];
rz(0.97719976) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7203622) q[1];
sx q[1];
rz(-2.1167397) q[1];
sx q[1];
rz(0.027688428) q[1];
x q[2];
rz(-0.67354789) q[3];
sx q[3];
rz(-1.827497) q[3];
sx q[3];
rz(1.7208769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6860883) q[2];
sx q[2];
rz(-2.4839165) q[2];
sx q[2];
rz(0.28798506) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.3214313) q[3];
sx q[3];
rz(0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31036723) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(-1.2121375) q[0];
rz(-0.37777004) q[1];
sx q[1];
rz(-1.2021844) q[1];
sx q[1];
rz(1.4935965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6389248) q[0];
sx q[0];
rz(-1.8687316) q[0];
sx q[0];
rz(-1.7940679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1475032) q[2];
sx q[2];
rz(-2.9967542) q[2];
sx q[2];
rz(2.7763979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32840604) q[1];
sx q[1];
rz(-1.9771736) q[1];
sx q[1];
rz(-2.3996668) q[1];
rz(1.3573523) q[3];
sx q[3];
rz(-0.5118013) q[3];
sx q[3];
rz(1.7970999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.2517694) q[2];
rz(-0.88636032) q[3];
sx q[3];
rz(-2.3754407) q[3];
sx q[3];
rz(-0.10929508) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-2.1003523) q[0];
sx q[0];
rz(-0.1782724) q[0];
rz(-3.0687304) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(1.1791621) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423791) q[0];
sx q[0];
rz(-1.8633435) q[0];
sx q[0];
rz(0.53036687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91782848) q[2];
sx q[2];
rz(-0.73790109) q[2];
sx q[2];
rz(0.81673056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0563911) q[1];
sx q[1];
rz(-1.8279652) q[1];
sx q[1];
rz(-2.609054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0696899) q[3];
sx q[3];
rz(-1.9074829) q[3];
sx q[3];
rz(-1.5130373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4420085) q[2];
sx q[2];
rz(-0.99094355) q[2];
sx q[2];
rz(-2.1098095) q[2];
rz(1.3423086) q[3];
sx q[3];
rz(-1.7248025) q[3];
sx q[3];
rz(2.9163196) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.306698) q[1];
sx q[1];
rz(-1.4790165) q[1];
sx q[1];
rz(-2.0526989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63554791) q[0];
sx q[0];
rz(-1.6354927) q[0];
sx q[0];
rz(-3.0788172) q[0];
rz(-pi) q[1];
rz(-0.18386545) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(1.136214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(1.0165434) q[1];
rz(-2.0652566) q[3];
sx q[3];
rz(-0.38443243) q[3];
sx q[3];
rz(-1.9264551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-2.4536295) q[2];
sx q[2];
rz(0.069996746) q[2];
rz(-2.2648515) q[3];
sx q[3];
rz(-1.8165959) q[3];
sx q[3];
rz(-0.35818067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.6246673) q[1];
sx q[1];
rz(-1.0472052) q[1];
sx q[1];
rz(-0.44890961) q[1];
rz(2.0786053) q[2];
sx q[2];
rz(-0.72401902) q[2];
sx q[2];
rz(0.25821092) q[2];
rz(-1.4528965) q[3];
sx q[3];
rz(-2.8040734) q[3];
sx q[3];
rz(0.87344195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
