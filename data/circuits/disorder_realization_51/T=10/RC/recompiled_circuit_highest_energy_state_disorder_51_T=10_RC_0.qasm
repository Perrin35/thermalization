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
rz(0.32938862) q[0];
sx q[0];
rz(3.7731054) q[0];
sx q[0];
rz(9.4256529) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54291081) q[0];
sx q[0];
rz(-3.0474159) q[0];
sx q[0];
rz(-1.3046632) q[0];
rz(-pi) q[1];
rz(1.4232741) q[2];
sx q[2];
rz(-1.9468687) q[2];
sx q[2];
rz(-3.0126115) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4507323) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(0.91207204) q[1];
x q[2];
rz(1.4963989) q[3];
sx q[3];
rz(-1.5870023) q[3];
sx q[3];
rz(0.61283535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22902809) q[2];
sx q[2];
rz(-1.8077069) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(-2.8424756) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107373) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(-0.77450007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82371432) q[0];
sx q[0];
rz(-1.5817989) q[0];
sx q[0];
rz(-1.3925793) q[0];
x q[1];
rz(2.108091) q[2];
sx q[2];
rz(-2.3024493) q[2];
sx q[2];
rz(-3.0947859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0424287) q[1];
sx q[1];
rz(-1.1981257) q[1];
sx q[1];
rz(1.8120239) q[1];
rz(-pi) q[2];
rz(-1.1011613) q[3];
sx q[3];
rz(-1.7723933) q[3];
sx q[3];
rz(-2.2045709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.896686) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(2.1265105) q[2];
rz(0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(-2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.86431137) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(-2.1306254) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(-3.0027622) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75051266) q[0];
sx q[0];
rz(-1.2475999) q[0];
sx q[0];
rz(-1.1884407) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9857668) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(2.6443554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0142747) q[1];
sx q[1];
rz(-1.7671314) q[1];
sx q[1];
rz(-1.8683968) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7458472) q[3];
sx q[3];
rz(-1.850046) q[3];
sx q[3];
rz(-1.6498914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6802754) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(-0.3581363) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-0.3048234) q[0];
rz(-0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-0.92794424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5261032) q[0];
sx q[0];
rz(-0.23357059) q[0];
sx q[0];
rz(2.156267) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60073845) q[2];
sx q[2];
rz(-1.033412) q[2];
sx q[2];
rz(-1.5058277) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.881262) q[1];
sx q[1];
rz(-1.4196383) q[1];
sx q[1];
rz(-0.69044729) q[1];
x q[2];
rz(2.4797012) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(-0.5058561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-0.18816571) q[2];
rz(-1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4417878) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(-3.1222043) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.2394989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0415619) q[0];
sx q[0];
rz(-1.8503555) q[0];
sx q[0];
rz(-2.1651405) q[0];
rz(-0.92186982) q[2];
sx q[2];
rz(-0.94266329) q[2];
sx q[2];
rz(-0.79337304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46249786) q[1];
sx q[1];
rz(-2.0789008) q[1];
sx q[1];
rz(1.2485882) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2438888) q[3];
sx q[3];
rz(-2.5335059) q[3];
sx q[3];
rz(1.9185818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(0.2462247) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(-2.4024409) q[0];
rz(-2.7958561) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(2.2703222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7428674) q[0];
sx q[0];
rz(-0.79116066) q[0];
sx q[0];
rz(1.6683116) q[0];
rz(-pi) q[1];
rz(2.5315782) q[2];
sx q[2];
rz(-0.84976174) q[2];
sx q[2];
rz(-0.80243669) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(1.9597998) q[1];
x q[2];
rz(2.0454117) q[3];
sx q[3];
rz(-1.0132917) q[3];
sx q[3];
rz(1.4589918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(2.2147307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0480014) q[0];
sx q[0];
rz(-2.8107852) q[0];
sx q[0];
rz(-0.76865102) q[0];
x q[1];
rz(0.92717391) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(1.3869029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95185223) q[1];
sx q[1];
rz(-1.2092672) q[1];
sx q[1];
rz(-0.23747634) q[1];
rz(-1.5038483) q[3];
sx q[3];
rz(-0.62426011) q[3];
sx q[3];
rz(-2.3926596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(-0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153862) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(-0.4183847) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(-0.13430886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90420049) q[0];
sx q[0];
rz(-2.2456894) q[0];
sx q[0];
rz(0.86752059) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1654794) q[2];
sx q[2];
rz(-1.7421054) q[2];
sx q[2];
rz(1.5660945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5843552) q[1];
sx q[1];
rz(-1.7533301) q[1];
sx q[1];
rz(1.4339158) q[1];
x q[2];
rz(2.0401434) q[3];
sx q[3];
rz(-1.3162287) q[3];
sx q[3];
rz(1.001844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(2.4995787) q[2];
rz(-1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5409656) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(1.8070096) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(-2.1122011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6717259) q[0];
sx q[0];
rz(-2.436536) q[0];
sx q[0];
rz(-1.6204349) q[0];
rz(-pi) q[1];
rz(-3.0156187) q[2];
sx q[2];
rz(-0.306189) q[2];
sx q[2];
rz(0.15737113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16949108) q[1];
sx q[1];
rz(-1.5268561) q[1];
sx q[1];
rz(1.8010587) q[1];
x q[2];
rz(0.179788) q[3];
sx q[3];
rz(-2.6670293) q[3];
sx q[3];
rz(-3.0347412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42334291) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(1.8839802) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(0.14258252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6194942) q[0];
sx q[0];
rz(-1.5987248) q[0];
sx q[0];
rz(-2.9832207) q[0];
rz(-1.8480186) q[2];
sx q[2];
rz(-2.1810594) q[2];
sx q[2];
rz(-2.0137613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0847454) q[1];
sx q[1];
rz(-1.0029625) q[1];
sx q[1];
rz(0.3155057) q[1];
rz(1.110465) q[3];
sx q[3];
rz(-1.8813731) q[3];
sx q[3];
rz(-0.47720695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3451781) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(0.26819116) q[2];
rz(-1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0161229) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(0.62025537) q[2];
sx q[2];
rz(-0.56369416) q[2];
sx q[2];
rz(-0.77947215) q[2];
rz(1.2760194) q[3];
sx q[3];
rz(-1.5073007) q[3];
sx q[3];
rz(-0.56137847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
