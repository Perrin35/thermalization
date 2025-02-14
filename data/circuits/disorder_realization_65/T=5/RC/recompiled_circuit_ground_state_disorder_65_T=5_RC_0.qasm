OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4545257) q[0];
sx q[0];
rz(-1.2626167) q[0];
sx q[0];
rz(0.73448056) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(1.8416815) q[1];
sx q[1];
rz(11.587974) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7534417) q[0];
sx q[0];
rz(-1.4343669) q[0];
sx q[0];
rz(3.1336954) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4327203) q[2];
sx q[2];
rz(-1.9908496) q[2];
sx q[2];
rz(1.5839029) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20767388) q[1];
sx q[1];
rz(-1.7130995) q[1];
sx q[1];
rz(1.2486904) q[1];
x q[2];
rz(-3.0618598) q[3];
sx q[3];
rz(-2.0679722) q[3];
sx q[3];
rz(-0.91547188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3413099) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(0.71551234) q[2];
rz(-1.4095151) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187123) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(3.1067644) q[0];
rz(-0.71827978) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(1.9504331) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5226591) q[0];
sx q[0];
rz(-2.1858523) q[0];
sx q[0];
rz(0.19490029) q[0];
x q[1];
rz(2.9600675) q[2];
sx q[2];
rz(-2.4511538) q[2];
sx q[2];
rz(1.4140364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0409567) q[1];
sx q[1];
rz(-1.2946048) q[1];
sx q[1];
rz(0.60224288) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19050772) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(-2.0634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8027773) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(3.0954933) q[2];
rz(-0.80373803) q[3];
sx q[3];
rz(-1.9019144) q[3];
sx q[3];
rz(0.55155915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4285202) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(0.61306104) q[0];
rz(-0.26519457) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(3.0934966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0408024) q[0];
sx q[0];
rz(-1.5682966) q[0];
sx q[0];
rz(-2.0526171) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5597495) q[2];
sx q[2];
rz(-2.0177236) q[2];
sx q[2];
rz(2.4012411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7386507) q[1];
sx q[1];
rz(-1.623282) q[1];
sx q[1];
rz(-1.6305292) q[1];
rz(1.223391) q[3];
sx q[3];
rz(-1.239353) q[3];
sx q[3];
rz(-3.0104266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6951822) q[2];
sx q[2];
rz(-1.9789275) q[2];
sx q[2];
rz(2.7282558) q[2];
rz(-0.6684331) q[3];
sx q[3];
rz(-1.8139402) q[3];
sx q[3];
rz(1.7641164) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17871019) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(2.9982153) q[0];
rz(-2.7168221) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(-0.43407789) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052267678) q[0];
sx q[0];
rz(-0.10073034) q[0];
sx q[0];
rz(1.1830551) q[0];
x q[1];
rz(2.4894478) q[2];
sx q[2];
rz(-1.6639478) q[2];
sx q[2];
rz(2.8217614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.886288) q[1];
sx q[1];
rz(-1.7105173) q[1];
sx q[1];
rz(2.0934859) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9346385) q[3];
sx q[3];
rz(-2.3194139) q[3];
sx q[3];
rz(1.5311629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8726337) q[2];
sx q[2];
rz(-0.7577529) q[2];
sx q[2];
rz(-0.41716519) q[2];
rz(-0.68552351) q[3];
sx q[3];
rz(-1.7020107) q[3];
sx q[3];
rz(-3.0912257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3454469) q[0];
sx q[0];
rz(-2.8247927) q[0];
sx q[0];
rz(1.5078804) q[0];
rz(-0.66868526) q[1];
sx q[1];
rz(-1.1983127) q[1];
sx q[1];
rz(-1.1315469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1567845) q[0];
sx q[0];
rz(-0.97301345) q[0];
sx q[0];
rz(3.1020107) q[0];
rz(0.3338694) q[2];
sx q[2];
rz(-1.0354932) q[2];
sx q[2];
rz(-1.4620505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3858429) q[1];
sx q[1];
rz(-1.7569225) q[1];
sx q[1];
rz(-0.47551861) q[1];
rz(-0.82238166) q[3];
sx q[3];
rz(-2.9677792) q[3];
sx q[3];
rz(-1.2863152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23182997) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(2.6798999) q[2];
rz(-2.8716904) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(-0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5676024) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(-2.0800166) q[0];
rz(2.4827982) q[1];
sx q[1];
rz(-1.4219204) q[1];
sx q[1];
rz(-0.40491358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.957307) q[0];
sx q[0];
rz(-0.71249639) q[0];
sx q[0];
rz(-0.16100968) q[0];
x q[1];
rz(1.986386) q[2];
sx q[2];
rz(-2.0743557) q[2];
sx q[2];
rz(-1.0114618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3337219) q[1];
sx q[1];
rz(-2.0757339) q[1];
sx q[1];
rz(-0.15926475) q[1];
rz(1.9642716) q[3];
sx q[3];
rz(-2.4654536) q[3];
sx q[3];
rz(2.2666988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.095526516) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(1.3859762) q[2];
rz(2.4447377) q[3];
sx q[3];
rz(-2.160852) q[3];
sx q[3];
rz(-2.3314886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54414576) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(-0.59462732) q[0];
rz(-1.1320629) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(2.6304257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0027255) q[0];
sx q[0];
rz(-2.1423233) q[0];
sx q[0];
rz(-1.6353428) q[0];
rz(-0.86894247) q[2];
sx q[2];
rz(-2.7232183) q[2];
sx q[2];
rz(-1.4792031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0446544) q[1];
sx q[1];
rz(-1.4739828) q[1];
sx q[1];
rz(-1.3083878) q[1];
x q[2];
rz(-2.0603544) q[3];
sx q[3];
rz(-2.2415906) q[3];
sx q[3];
rz(1.2257934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96659294) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(0.28755507) q[2];
rz(-1.1047085) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(-3.0159359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6519046) q[0];
sx q[0];
rz(-2.7629485) q[0];
sx q[0];
rz(0.65377671) q[0];
rz(-2.6405624) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(-2.3562866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619072) q[0];
sx q[0];
rz(-2.5568049) q[0];
sx q[0];
rz(0.25081046) q[0];
x q[1];
rz(-2.5748524) q[2];
sx q[2];
rz(-1.9575685) q[2];
sx q[2];
rz(-2.6865328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4738087) q[1];
sx q[1];
rz(-0.89783339) q[1];
sx q[1];
rz(-1.8401309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9124746) q[3];
sx q[3];
rz(-1.6506628) q[3];
sx q[3];
rz(0.95019482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.952897) q[2];
sx q[2];
rz(-1.6152103) q[2];
sx q[2];
rz(0.3696028) q[2];
rz(0.26436198) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0109176) q[0];
sx q[0];
rz(-0.70931804) q[0];
sx q[0];
rz(1.6226907) q[0];
rz(-1.9021289) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(-2.1942203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560871) q[0];
sx q[0];
rz(-1.5099021) q[0];
sx q[0];
rz(1.1709309) q[0];
rz(1.6336254) q[2];
sx q[2];
rz(-1.9312973) q[2];
sx q[2];
rz(-2.6897827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92129642) q[1];
sx q[1];
rz(-1.7361169) q[1];
sx q[1];
rz(-2.9784747) q[1];
rz(-1.2227433) q[3];
sx q[3];
rz(-2.9677941) q[3];
sx q[3];
rz(2.9617975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6685247) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(-1.1119941) q[2];
rz(-2.9396577) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(-0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875459) q[0];
sx q[0];
rz(-0.18809479) q[0];
sx q[0];
rz(-1.6756219) q[0];
rz(1.5415883) q[1];
sx q[1];
rz(-1.8406638) q[1];
sx q[1];
rz(0.25513729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2422243) q[0];
sx q[0];
rz(-1.5535019) q[0];
sx q[0];
rz(1.4370741) q[0];
rz(-pi) q[1];
rz(1.8996283) q[2];
sx q[2];
rz(-0.82482289) q[2];
sx q[2];
rz(-2.0212295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0560811) q[1];
sx q[1];
rz(-1.8485118) q[1];
sx q[1];
rz(3.059832) q[1];
x q[2];
rz(0.76948036) q[3];
sx q[3];
rz(-1.7090685) q[3];
sx q[3];
rz(-1.4797803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5280891) q[2];
sx q[2];
rz(-1.9518423) q[2];
sx q[2];
rz(2.7733754) q[2];
rz(-2.8585989) q[3];
sx q[3];
rz(-0.26809433) q[3];
sx q[3];
rz(-2.8158269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1267283) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(-2.6425843) q[1];
sx q[1];
rz(-1.2393163) q[1];
sx q[1];
rz(-2.7543482) q[1];
rz(0.3327315) q[2];
sx q[2];
rz(-1.4495779) q[2];
sx q[2];
rz(-2.4548893) q[2];
rz(1.7718689) q[3];
sx q[3];
rz(-1.5101505) q[3];
sx q[3];
rz(2.8240614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
