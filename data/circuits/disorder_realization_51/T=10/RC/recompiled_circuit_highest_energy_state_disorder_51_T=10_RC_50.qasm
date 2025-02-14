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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76287718) q[0];
sx q[0];
rz(-1.595531) q[0];
sx q[0];
rz(-1.6616761) q[0];
x q[1];
rz(2.7617747) q[2];
sx q[2];
rz(-1.7079412) q[2];
sx q[2];
rz(1.4963381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7172437) q[1];
sx q[1];
rz(-2.1645438) q[1];
sx q[1];
rz(-0.51148606) q[1];
rz(-pi) q[2];
rz(-1.3561068) q[3];
sx q[3];
rz(-3.0654538) q[3];
sx q[3];
rz(-2.3977181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9125646) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(-0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107373) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(-0.99288565) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(0.77450007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555154) q[0];
sx q[0];
rz(-0.17855274) q[0];
sx q[0];
rz(-1.6327842) q[0];
rz(2.6235145) q[2];
sx q[2];
rz(-0.87730125) q[2];
sx q[2];
rz(2.3665646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0424287) q[1];
sx q[1];
rz(-1.1981257) q[1];
sx q[1];
rz(1.3295688) q[1];
rz(-pi) q[2];
rz(-2.9162972) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(2.6090906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2449067) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(2.1306254) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(-0.13883042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9897744) q[0];
sx q[0];
rz(-0.4954557) q[0];
sx q[0];
rz(0.83924967) q[0];
rz(-2.6723318) q[2];
sx q[2];
rz(-0.77232328) q[2];
sx q[2];
rz(-2.0283716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0142747) q[1];
sx q[1];
rz(-1.3744613) q[1];
sx q[1];
rz(1.2731958) q[1];
x q[2];
rz(2.5957727) q[3];
sx q[3];
rz(-2.8132317) q[3];
sx q[3];
rz(2.2203746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(0.67874587) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-1.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(0.3048234) q[0];
rz(-2.7335956) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-2.2136484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5281677) q[0];
sx q[0];
rz(-1.6990464) q[0];
sx q[0];
rz(-1.7665461) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60073845) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(1.5058277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18670652) q[1];
sx q[1];
rz(-0.88973239) q[1];
sx q[1];
rz(-1.3757371) q[1];
rz(-pi) q[2];
rz(-0.66189142) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(-0.5058561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(0.18816571) q[2];
rz(1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(1.060932) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.2394989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4864522) q[0];
sx q[0];
rz(-2.1391272) q[0];
sx q[0];
rz(0.33354946) q[0];
rz(2.4024348) q[2];
sx q[2];
rz(-2.0817167) q[2];
sx q[2];
rz(-0.3581274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0026999) q[1];
sx q[1];
rz(-0.59400105) q[1];
sx q[1];
rz(0.51704375) q[1];
rz(-pi) q[2];
rz(-0.89770384) q[3];
sx q[3];
rz(-0.60808676) q[3];
sx q[3];
rz(-1.9185818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1218607) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(-1.8149553) q[2];
rz(-0.2462247) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-2.2897913) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(0.34573653) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5369878) q[0];
sx q[0];
rz(-2.3571627) q[0];
sx q[0];
rz(3.0434199) q[0];
x q[1];
rz(-2.3910693) q[2];
sx q[2];
rz(-1.125968) q[2];
sx q[2];
rz(1.9409279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.516219) q[1];
sx q[1];
rz(-1.9587079) q[1];
sx q[1];
rz(-0.078887786) q[1];
rz(2.5302251) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(-0.37722019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(0.87257067) q[2];
rz(-1.568659) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0503814) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(-0.92686191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2970718) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(-1.336444) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71452801) q[2];
sx q[2];
rz(-2.0863669) q[2];
sx q[2];
rz(2.553568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7898579) q[1];
sx q[1];
rz(-2.7119293) q[1];
sx q[1];
rz(-2.1273145) q[1];
rz(0.94759946) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(-2.2653803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0285792) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(-2.8431559) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153862) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(-2.723208) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1100562) q[0];
sx q[0];
rz(-0.93288619) q[0];
sx q[0];
rz(-0.67968926) q[0];
x q[1];
rz(-2.9555126) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(-3.0732791) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5843552) q[1];
sx q[1];
rz(-1.3882625) q[1];
sx q[1];
rz(-1.7076769) q[1];
rz(-pi) q[2];
rz(-0.28388309) q[3];
sx q[3];
rz(-1.1177269) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-1.0293915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6717259) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(-1.6204349) q[0];
rz(2.8376828) q[2];
sx q[2];
rz(-1.5329156) q[2];
sx q[2];
rz(1.5335976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9255213) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(-1.3804803) q[1];
rz(-1.6624032) q[3];
sx q[3];
rz(-2.0370968) q[3];
sx q[3];
rz(-2.8331851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42334291) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(-0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-2.7227962) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(0.14258252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0973546) q[0];
sx q[0];
rz(-1.4124866) q[0];
sx q[0];
rz(1.5990785) q[0];
x q[1];
rz(0.62866399) q[2];
sx q[2];
rz(-1.3445879) q[2];
sx q[2];
rz(0.60461254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0847454) q[1];
sx q[1];
rz(-1.0029625) q[1];
sx q[1];
rz(-2.826087) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94513388) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(-1.6459203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0161229) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(1.9228946) q[2];
sx q[2];
rz(-2.0205971) q[2];
sx q[2];
rz(-1.481075) q[2];
rz(-1.3553452) q[3];
sx q[3];
rz(-2.8402495) q[3];
sx q[3];
rz(0.80339669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
