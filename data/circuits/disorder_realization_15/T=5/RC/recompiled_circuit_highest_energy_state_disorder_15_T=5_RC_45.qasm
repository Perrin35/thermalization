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
rz(1.6871356) q[0];
sx q[0];
rz(1.9782826) q[0];
sx q[0];
rz(10.692698) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(2.9719404) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21714144) q[0];
sx q[0];
rz(-1.0944774) q[0];
sx q[0];
rz(-2.7360271) q[0];
x q[1];
rz(1.1450139) q[2];
sx q[2];
rz(-2.5822431) q[2];
sx q[2];
rz(2.1970791) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4028266) q[1];
sx q[1];
rz(-1.9662274) q[1];
sx q[1];
rz(0.02301245) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98244169) q[3];
sx q[3];
rz(-0.37068493) q[3];
sx q[3];
rz(1.4511758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9156645) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(2.9018371) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(-0.83893004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0073485) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(2.0452621) q[0];
rz(1.6657375) q[1];
sx q[1];
rz(-1.219039) q[1];
sx q[1];
rz(-2.347167) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61928328) q[0];
sx q[0];
rz(-1.2566152) q[0];
sx q[0];
rz(-2.7854325) q[0];
rz(-pi) q[1];
rz(-2.3874823) q[2];
sx q[2];
rz(-2.0144785) q[2];
sx q[2];
rz(-0.59690969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5880044) q[1];
sx q[1];
rz(-2.70397) q[1];
sx q[1];
rz(-0.2511843) q[1];
rz(-pi) q[2];
rz(1.3737455) q[3];
sx q[3];
rz(-2.3757977) q[3];
sx q[3];
rz(-1.9869508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0558001) q[2];
sx q[2];
rz(-1.5689359) q[2];
sx q[2];
rz(-2.9158578) q[2];
rz(0.24383946) q[3];
sx q[3];
rz(-2.1832681) q[3];
sx q[3];
rz(-2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589515) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(0.42695811) q[0];
rz(-2.7985165) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-2.8899736) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8247525) q[0];
sx q[0];
rz(-1.7946294) q[0];
sx q[0];
rz(0.59455183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22253665) q[2];
sx q[2];
rz(-2.3250569) q[2];
sx q[2];
rz(2.726462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7921264) q[1];
sx q[1];
rz(-1.9053962) q[1];
sx q[1];
rz(-1.7484596) q[1];
rz(-0.050046845) q[3];
sx q[3];
rz(-2.280683) q[3];
sx q[3];
rz(-2.0662226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.130927) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(-2.6962386) q[2];
rz(1.917786) q[3];
sx q[3];
rz(-1.2659975) q[3];
sx q[3];
rz(-2.7538917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7447516) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(-1.9388439) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(-2.4045827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4241959) q[0];
sx q[0];
rz(-1.5079466) q[0];
sx q[0];
rz(-0.24817384) q[0];
rz(0.36815181) q[2];
sx q[2];
rz(-2.2128004) q[2];
sx q[2];
rz(-2.8411691) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1070246) q[1];
sx q[1];
rz(-0.66063297) q[1];
sx q[1];
rz(-0.28566499) q[1];
x q[2];
rz(2.9748671) q[3];
sx q[3];
rz(-1.5136216) q[3];
sx q[3];
rz(2.056332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0143339) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(0.8117525) q[2];
rz(0.74770606) q[3];
sx q[3];
rz(-2.2816608) q[3];
sx q[3];
rz(-0.060613304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6478445) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(1.3519979) q[0];
rz(-2.3494675) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(-2.1314714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57841741) q[0];
sx q[0];
rz(-2.0978598) q[0];
sx q[0];
rz(-0.38946797) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5158122) q[2];
sx q[2];
rz(-1.1786945) q[2];
sx q[2];
rz(-0.34811172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0466442) q[1];
sx q[1];
rz(-1.4913591) q[1];
sx q[1];
rz(0.55752505) q[1];
rz(1.406448) q[3];
sx q[3];
rz(-1.6858628) q[3];
sx q[3];
rz(3.1361406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1048364) q[2];
sx q[2];
rz(-2.3827621) q[2];
sx q[2];
rz(-1.7581615) q[2];
rz(-0.029959921) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(-1.8647319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5112011) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(0.1314441) q[0];
rz(-0.99880544) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.7703895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4192062) q[0];
sx q[0];
rz(-0.50438577) q[0];
sx q[0];
rz(2.3729411) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7741469) q[2];
sx q[2];
rz(-2.2468745) q[2];
sx q[2];
rz(-2.072056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2590717) q[1];
sx q[1];
rz(-1.4132199) q[1];
sx q[1];
rz(-1.8474008) q[1];
rz(-pi) q[2];
rz(-0.71089069) q[3];
sx q[3];
rz(-1.5681453) q[3];
sx q[3];
rz(-0.34429911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78352952) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-1.0895458) q[2];
rz(2.774636) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628196) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(-1.577042) q[1];
sx q[1];
rz(-1.4955474) q[1];
sx q[1];
rz(0.73211342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0769147) q[0];
sx q[0];
rz(-0.87192649) q[0];
sx q[0];
rz(-1.3069606) q[0];
rz(-pi) q[1];
rz(0.15190923) q[2];
sx q[2];
rz(-2.4503802) q[2];
sx q[2];
rz(2.76413) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.435675) q[1];
sx q[1];
rz(-0.49723782) q[1];
sx q[1];
rz(1.2364718) q[1];
rz(1.7978396) q[3];
sx q[3];
rz(-1.4770835) q[3];
sx q[3];
rz(0.39484398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6096036) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(0.58094376) q[2];
rz(-1.1254492) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(-0.17053764) q[0];
rz(-1.2455617) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(0.68444288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67407284) q[0];
sx q[0];
rz(-2.0430123) q[0];
sx q[0];
rz(-1.945482) q[0];
rz(-pi) q[1];
rz(-0.48284097) q[2];
sx q[2];
rz(-1.6833268) q[2];
sx q[2];
rz(-1.4275328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43757968) q[1];
sx q[1];
rz(-1.8270758) q[1];
sx q[1];
rz(-1.8248454) q[1];
rz(-pi) q[2];
rz(-2.4926458) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(-3.0779882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.838852) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(-1.1203241) q[2];
rz(2.4773347) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(-0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(0.63825178) q[0];
rz(3.1328746) q[1];
sx q[1];
rz(-2.0254841) q[1];
sx q[1];
rz(-2.7353824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531949) q[0];
sx q[0];
rz(-2.8450845) q[0];
sx q[0];
rz(1.1738846) q[0];
rz(-pi) q[1];
rz(-1.4582602) q[2];
sx q[2];
rz(-0.43607682) q[2];
sx q[2];
rz(1.2051518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49055028) q[1];
sx q[1];
rz(-2.3458651) q[1];
sx q[1];
rz(2.9814238) q[1];
rz(-pi) q[2];
rz(-3.0681455) q[3];
sx q[3];
rz(-1.5454419) q[3];
sx q[3];
rz(3.1318093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.093988769) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(0.36063933) q[2];
rz(-2.4141198) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(-0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(2.9720921) q[0];
rz(0.70973474) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(1.08606) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5026131) q[0];
sx q[0];
rz(-1.3492246) q[0];
sx q[0];
rz(0.53934877) q[0];
x q[1];
rz(0.96699826) q[2];
sx q[2];
rz(-2.6259086) q[2];
sx q[2];
rz(0.57920757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0154961) q[1];
sx q[1];
rz(-1.0119146) q[1];
sx q[1];
rz(-0.11591594) q[1];
rz(-2.7602155) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(-1.7538296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3544932) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(-1.1981111) q[3];
sx q[3];
rz(-1.5912278) q[3];
sx q[3];
rz(-0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030466) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(-1.5003149) q[2];
sx q[2];
rz(-1.6432091) q[2];
sx q[2];
rz(2.5753449) q[2];
rz(1.4403338) q[3];
sx q[3];
rz(-1.4080019) q[3];
sx q[3];
rz(-0.27346591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
